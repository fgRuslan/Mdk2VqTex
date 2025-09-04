#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <random>

#include "decompressor.h"

#pragma pack(push,1)
struct TextureHeader {
    uint32_t magic = 0x000007D1;
    uint32_t version = 0x00000015;
    uint32_t width;
    uint32_t height;
    uint32_t flags = 0;
    uint32_t data_size;
};
#pragma pack(pop)

// --- Lookup Tables (Values from original code) ---
const uint32_t dword_4B3B80[] = { 0x40, 0x5E };
const uint32_t dword_4B3B88[] = { 0, 0x20 };
const uint32_t dword_4B3B90[] = { 0, 0x20 };

// Utility structures
struct Color {
    union {
        struct {
            uint8_t r, g, b, a;
        };
        uint8_t components[4];
    };
	
	uint8_t operator[](int i)
	{
		return this->components[i];
	}
};

struct ColorF {
    float r, g, b, a;
};

// Precompute quantization tables for 5-bit to 8-bit and reverse
uint8_t quant5[256];
uint8_t expand5[32];
void init_quant_tables() {
    for (int q = 0; q < 32; ++q) {
        expand5[q] = (q << 3) | (q >> 2);
    }
    for (int v = 0; v < 256; ++v) {
        int best_q = 0;
        int best_err = std::numeric_limits<int>::max();
        for (int q = 0; q < 32; ++q) {
            int err = std::abs(static_cast<int>(expand5[q]) - v);
            if (err < best_err) {
                best_err = err;
                best_q = q;
            }
        }
        quant5[v] = best_q;
    }
}

// Bit setter
inline void set_bit(std::vector<uint32_t>& data, int bit_index, uint32_t value) {
    if ((bit_index >> 5) >= data.size()) return;
    if (value) {
        data[bit_index >> 5] |= (1u << (bit_index & 0x1F));
    } else {
        data[bit_index >> 5] &= ~(1u << (bit_index & 0x1F));
    }
}

// Distance between colors
float color_dist(const Color& c1, const Color& c2, bool include_alpha) {
    float dr = c1.r - c2.r;
    float dg = c1.g - c2.g;
    float db = c1.b - c2.b;
    float da = include_alpha ? c1.a - c2.a : 0;
    return dr*dr + dg*dg + db*db + da*da;
}

// Extract 4x8 block from image
void extract_block(const unsigned char* image, int width, int height, int block_x, int block_y, Color block_pixels[4][8]) {
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 8; ++col) {
            int px = block_x * 8 + col;
            int py = block_y * 4 + row;
            if (px < width && py < height) {
                int offset = (py * width + px) * 4;
                // Swap R and B to match internal representation
                block_pixels[row][col] = {image[offset + 2], image[offset + 1], image[offset + 0], image[offset + 3]};
            } else {
                block_pixels[row][col] = {0, 0, 0, 0};
            }
        }
    }
}

// Compute SSE error by decompressing and comparing
float compute_block_error(Color original[4][8], const std::vector<uint32_t>& block_data) {
    uint32_t temp_buffer[4 * 8] = {0};
    decode_block(8, 0, 0, temp_buffer, block_data); // Use width=8 for small buffer
    float sse = 0.0f;
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 8; ++col) {
            uint32_t p = temp_buffer[row * 8 + col];
            Color recon = {static_cast<uint8_t>(p & 0xFF), static_cast<uint8_t>((p >> 8) & 0xFF),
                          static_cast<uint8_t>((p >> 16) & 0xFF), static_cast<uint8_t>((p >> 24) & 0xFF)};
            sse += color_dist(original[row][col], recon, true);
        }
    }
    return sse;
}

// Find max distance endpoints
void find_endpoints(const std::vector<Color>& pixels, Color& end1, Color& end2, bool include_alpha) {
    float max_dist = -1.0f;
    for (const auto& p1 : pixels) {
        for (const auto& p2 : pixels) {
            float d = color_dist(p1, p2, include_alpha);
            if (d > max_dist) {
                max_dist = d;
                end1 = p1;
                end2 = p2;
            }
        }
    }
}

// Simple k-means for clustering
void kmeans_cluster(const std::vector<Color>& pixels, int k, Color centers[4], bool include_alpha) {
    if (pixels.empty()) return;
    std::vector<ColorF> cf_centers(k);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, pixels.size() - 1);
    for (int i = 0; i < k; ++i) {
        Color p = pixels[dis(gen)];
        cf_centers[i] = {static_cast<float>(p.r), static_cast<float>(p.g), static_cast<float>(p.b), static_cast<float>(p.a)};
    }
    for (int iter = 0; iter < 10; ++iter) {
        std::vector<std::vector<ColorF>> groups(k);
        for (const auto& p : pixels) {
            ColorF pf = {static_cast<float>(p.r), static_cast<float>(p.g), static_cast<float>(p.b), static_cast<float>(p.a)};
            int best = 0;
            Color c0 = {static_cast<uint8_t>(cf_centers[0].r), static_cast<uint8_t>(cf_centers[0].g), static_cast<uint8_t>(cf_centers[0].b), static_cast<uint8_t>(cf_centers[0].a)};
            float best_d = color_dist(p, c0, include_alpha);
            for (int j = 1; j < k; ++j) {
                Color cj = {static_cast<uint8_t>(cf_centers[j].r), static_cast<uint8_t>(cf_centers[j].g), static_cast<uint8_t>(cf_centers[j].b), static_cast<uint8_t>(cf_centers[j].a)};
                float d = color_dist(p, cj, include_alpha);
                if (d < best_d) {
                    best_d = d;
                    best = j;
                }
            }
            groups[best].push_back(pf);
        }
        for (int j = 0; j < k; ++j) {
            if (groups[j].empty()) continue;
            ColorF sum = {0,0,0,0};
            for (const auto& s : groups[j]) {
                sum.r += s.r;
                sum.g += s.g;
                sum.b += s.b;
                sum.a += s.a;
            }
            int sz = groups[j].size();
            sum.r /= sz; sum.g /= sz; sum.b /= sz; sum.a /= sz;
            cf_centers[j] = sum;
        }
    }
    for (int j = 0; j < k; ++j) {
        centers[j].components[0] = expand5[quant5[static_cast<uint8_t>(std::round(cf_centers[j].r))]];
        centers[j].components[1] = expand5[quant5[static_cast<uint8_t>(std::round(cf_centers[j].g))]];
        centers[j].components[2] = expand5[quant5[static_cast<uint8_t>(std::round(cf_centers[j].b))]];
        centers[j].components[3] = expand5[quant5[static_cast<uint8_t>(std::round(cf_centers[j].a))]];
    }
}

// Compress a single block
std::vector<uint32_t> compress_block(Color block_pixels[4][8]) {
    std::vector<uint32_t> best_data(4, 0);
    float best_error = std::numeric_limits<float>::max();

    bool has_transparency = false;
    std::vector<Color> all_pixels;
    std::vector<Color> opaque_pixels;
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 8; ++col) {
            all_pixels.push_back(block_pixels[row][col]);
            if (block_pixels[row][col].a < 128) {
                has_transparency = true;
            } else {
                opaque_pixels.push_back(block_pixels[row][col]);
            }
        }
    }

    // --- OPAQUE-ONLY MODES ---
    // These modes are only tested if the block is fully opaque.
    if (!has_transparency) {
        // --- Try Mode A (Opaque Variant, v27=0) ---
        {
            std::vector<uint32_t> data(4, 0);
            uint32_t a8 = 0x80000000; // v27=0
            
            for (int block_idx = 0; block_idx < 2; ++block_idx) {
                std::vector<Color> sub_pixels;
                for (int r = 0; r < 4; ++r) for (int c = 0; c < 4; ++c) sub_pixels.push_back(block_pixels[r][block_idx * 4 + c]);
                
                Color c1, c2;
                find_endpoints(sub_pixels, c1, c2, false);
                c1.r = expand5[quant5[c1.r]]; c1.g = expand5[quant5[c1.g]]; c1.b = expand5[quant5[c1.b]]; c1.a = 255;
                c2.r = expand5[quant5[c2.r]]; c2.g = expand5[quant5[c2.g]]; c2.b = expand5[quant5[c2.b]]; c2.a = 255;

                uint32_t color_base_offset = dword_4B3B80[block_idx];
                for (int i=0; i<4; ++i) {
                    uint32_t q1 = quant5[c1.components[i]], q2 = quant5[c2.components[i]];
                    for (int bit=0; bit<5; ++bit) {
                        set_bit(data, color_base_offset + i*5 + bit, (q1 >> bit) & 1);
                        set_bit(data, color_base_offset + 15 + i*5 + bit, (q2 >> bit) & 1);
                    }
                }
                if (block_idx == 0) a8 |= (0 << 29); else a8 |= (0 << 30);
                set_bit(data, dword_4B3B88[block_idx] + 1, 0);

                uint32_t palette[16];
                uint32_t c1_comps[4]={c1.r,c1.g,c1.b,c1.a}, c2_comps[4]={c2.r,c2.g,c2.b,c2.a};
                generate_palette_from_ida(palette, c1_comps, c2_comps, false);
                
                int bit_stream_offset = dword_4B3B88[block_idx];
                for (int r=0; r<4; ++r) for (int c=0; c<4; ++c) {
                    Color p = sub_pixels[r*4+c];
                    int best_idx = 0;
                    float min_d = std::numeric_limits<float>::max();
                    for (int i=0; i<4; ++i) {
                        Color pal_c = {(uint8_t)palette[i*4],(uint8_t)palette[i*4+1],(uint8_t)palette[i*4+2],(uint8_t)palette[i*4+3]};
                        float d = color_dist(p, pal_c, false);
                        if (d < min_d) { min_d = d; best_idx = i; }
                    }
                    int bit_offset = r * 8 + c * 2;
                    set_bit(data, bit_stream_offset + bit_offset, best_idx & 1);
                    set_bit(data, bit_stream_offset + bit_offset + 1, (best_idx >> 1) & 1);
                }
            }
            data[3] = a8;
            float err = compute_block_error(block_pixels, data);
            if (err < best_error) { best_error = err; best_data = data; }
        }

        // --- Try Mode C ---
        {
            std::vector<uint32_t> data(4, 0);
            data[3] = 0x40000000;
            Color centers[4];
            kmeans_cluster(opaque_pixels, 4, centers, false);
            for (int comp=0; comp<3; ++comp) {
                for (int bit=0; bit<5; ++bit) {
                    set_bit(data, 79+comp*5-15+bit, (quant5[centers[0].components[comp]]>>bit)&1);
                    set_bit(data, 79+comp*5+bit, (quant5[centers[1].components[comp]]>>bit)&1);
                    set_bit(data, 79+comp*5+15+bit, (quant5[centers[2].components[comp]]>>bit)&1);
                    set_bit(data, 79+comp*5+30+bit, (quant5[centers[3].components[comp]]>>bit)&1);
                }
            }
            for (int r=0; r<4; ++r) for (int c=0; c<8; ++c) {
                Color p = block_pixels[r][c];
                int best_idx = 0;
                float min_d = std::numeric_limits<float>::max();
                for (int i=0; i<4; ++i) {
                    float d = color_dist(p, centers[i], false);
                    if (d < min_d) { min_d = d; best_idx = i; }
                }
                int bit_offset = (c < 4 ? 0 : 0x20) + 2 * (r * 4 + (c & 3));
                set_bit(data, bit_offset, best_idx & 1);
                set_bit(data, bit_offset+1, (best_idx >> 1) & 1);
            }
            float err = compute_block_error(block_pixels, data);
            if (err < best_error) { best_error = err; best_data = data; }
        }
    }

    // --- ALPHA-SUPPORTING MODES ---
    // These modes are always tested.

    // --- Try Mode A (Transparent Variant, v27=1) ---
    {
        std::vector<uint32_t> data(4, 0);
        uint32_t a8 = 0x80000000 | (1 << 28); // v27=1
        
        for (int block_idx = 0; block_idx < 2; ++block_idx) {
            std::vector<Color> sub_pixels;
            std::vector<Color> opaque_sub_pixels;
            for (int r = 0; r < 4; ++r) for (int c = 0; c < 4; ++c) {
                Color p = block_pixels[r][block_idx * 4 + c];
                sub_pixels.push_back(p);
                if (p.a >= 128) opaque_sub_pixels.push_back(p);
            }
            
            Color c1, c2;
            if (opaque_sub_pixels.empty()) { // Safe fallback for fully transparent sub-blocks
                c1 = {255,255,255,255}; c2 = {0,0,0,255};
            } else {
                find_endpoints(opaque_sub_pixels, c1, c2, false);
            }

            c1.r = expand5[quant5[c1.r]]; c1.g = expand5[quant5[c1.g]]; c1.b = expand5[quant5[c1.b]]; c1.a = 255;
            c2.r = expand5[quant5[c2.r]]; c2.g = expand5[quant5[c2.g]]; c2.b = expand5[quant5[c2.b]]; c2.a = 255;

            uint32_t color_base_offset = dword_4B3B80[block_idx];
            for (int i=0; i<4; ++i) {
                uint32_t q1 = quant5[c1.components[i]], q2 = quant5[c2.components[i]];
                for (int bit=0; bit<5; ++bit) {
                    set_bit(data, color_base_offset + i*5 + bit, (q1 >> bit) & 1);
                    set_bit(data, color_base_offset + 15 + i*5 + bit, (q2 >> bit) & 1);
                }
            }
            if (block_idx == 0) a8 |= (0 << 29); else a8 |= (0 << 30);

            uint32_t palette[16];
            uint32_t c1_comps[4]={c1.r,c1.g,c1.b,c1.a}, c2_comps[4]={c2.r,c2.g,c2.b,c2.a};
            generate_palette_from_ida(palette, c1_comps, c2_comps, true);
            
            int bit_stream_offset = dword_4B3B88[block_idx];
            for (int r=0; r<4; ++r) for (int c=0; c<4; ++c) {
                Color p = sub_pixels[r*4+c];
                int best_idx = 0;
                if (p.a < 128) {
                    best_idx = 3;
                } else {
                    float min_d = std::numeric_limits<float>::max();
                    for (int i=0; i<3; ++i) { // Only check the 3 opaque colors
                        Color pal_c = {(uint8_t)palette[i*4],(uint8_t)palette[i*4+1],(uint8_t)palette[i*4+2],(uint8_t)palette[i*4+3]};
                        float d = color_dist(p, pal_c, false);
                        if (d < min_d) { min_d = d; best_idx = i; }
                    }
                }
                int bit_offset = r * 8 + c * 2;
                set_bit(data, bit_stream_offset + bit_offset, best_idx & 1);
                set_bit(data, bit_stream_offset + bit_offset + 1, (best_idx >> 1) & 1);
            }
        }
        data[3] = a8;
        float err = compute_block_error(block_pixels, data);
        if (err < best_error) { best_error = err; best_data = data; }
    }

    // --- Try Mode B ---
    {
        std::vector<uint32_t> data(4, 0);
        data[3] = 0;
        Color c1, c2;
        if (opaque_pixels.empty()) { // Safe fallback
            c1 = {255,255,255,255}; c2 = {0,0,0,255};
        } else {
            find_endpoints(opaque_pixels, c1, c2, false);
        }
        c1.r = expand5[quant5[c1.r]]; c1.g = expand5[quant5[c1.g]]; c1.b = expand5[quant5[c1.b]];
        c2.r = expand5[quant5[c2.r]]; c2.g = expand5[quant5[c2.g]]; c2.b = expand5[quant5[c2.b]];
        
        for (int i=0; i<3; ++i) {
            uint32_t q1 = quant5[c1.components[i]], q2 = quant5[c2.components[i]];
            for (int bit=0; bit<5; ++bit) {
                set_bit(data, 96 + i*5 + bit, (q1 >> bit) & 1);
                set_bit(data, 111 + i*5 + bit, (q2 >> bit) & 1);
            }
        }

        uint32_t palette[32];
        uint32_t c1_comps_rgb[3]={c1.r,c1.g,c1.b}, c2_comps_rgb[3]={c2.r,c2.g,c2.b};
        generate_palette_mode1(palette, c1_comps_rgb, c2_comps_rgb);

        for (int r=0; r<4; ++r) for (int c=0; c<8; ++c) {
            Color p = block_pixels[r][c];
            int best_idx = 0;
            if (p.a < 128) {
                best_idx = 7;
            } else {
                float min_d = std::numeric_limits<float>::max();
                for (int i=0; i<7; ++i) {
                    Color pal_c = {(uint8_t)palette[i*4],(uint8_t)palette[i*4+1],(uint8_t)palette[i*4+2],(uint8_t)palette[i*4+3]};
                    float d = color_dist(p, pal_c, false);
                    if (d < min_d) { min_d = d; best_idx = i; }
                }
            }
            int bit_offset = (c < 4 ? 0 : 0x30) + 3 * (r * 4 + (c & 3));
            set_bit(data, bit_offset, best_idx & 1);
            set_bit(data, bit_offset+1, (best_idx >> 1) & 1);
            set_bit(data, bit_offset+2, (best_idx >> 2) & 1);
        }
        float err = compute_block_error(block_pixels, data);
        if (err < best_error) { best_error = err; best_data = data; }
    }

    // --- Try Mode D (True Alpha Mode) ---
    for (int v58 = 0; v58 < 2; ++v58) {
        if (v58 == 0 && !has_transparency) continue;

        std::vector<uint32_t> data(4, 0);
        uint32_t a8 = 0x60000000 | (static_cast<uint32_t>(v58) << 28);
        
        Color c[3];
        if (v58 == 0) { // Explicit RGBA colors
            kmeans_cluster(all_pixels, 3, c, true);
        } else { // Interpolated RGBA colors
            std::vector<Color> left_p, right_p;
            for(int r=0; r<4; ++r) for(int col=0; col<4; ++col) left_p.push_back(block_pixels[r][col]);
            for(int r=0; r<4; ++r) for(int col=0; col<4; ++col) right_p.push_back(block_pixels[r][col+4]);
            Color l0, l1, r0, r1;
            if(left_p.empty()){ l0={0,0,0,0}; l1={0,0,0,0}; } else { find_endpoints(left_p, l0, l1, true); }
            if(right_p.empty()){ r0={0,0,0,0}; r1={0,0,0,0}; } else { find_endpoints(right_p, r0, r1, true); }
            c[0] = l0; c[2] = r0;
            c[1].r = (l1.r + r1.r) / 2; c[1].g = (l1.g + r1.g) / 2;
            c[1].b = (l1.b + r1.b) / 2; c[1].a = (l1.a + r1.a) / 2;
        }

        for(int i=0; i<3; ++i) for(int j=0; j<4; ++j) c[i].components[j] = expand5[quant5[c[i].components[j]]];

        for (int comp=0; comp<3; ++comp) for (int bit=0; bit<5; ++bit) {
            set_bit(data, 79+comp*5-15+bit, (quant5[c[0].components[comp]]>>bit)&1);
            set_bit(data, 79+comp*5+bit, (quant5[c[1].components[comp]]>>bit)&1);
            set_bit(data, 79+comp*5+15+bit, (quant5[c[2].components[comp]]>>bit)&1);
        }
        for (int bit=0; bit<5; ++bit) {
            set_bit(data, 96+13+bit, (quant5[c[0].a]>>bit)&1);
            set_bit(data, 101+13+bit, (quant5[c[1].a]>>bit)&1);
            set_bit(data, 106+13+bit, (quant5[c[2].a]>>bit)&1);
        }

        if(v58 == 0) {
            for (int r=0; r<4; ++r) for (int col=0; col<8; ++col) {
                Color p = block_pixels[r][col];
                int best_idx = 3;
                if (p.a >= 128) {
                    float min_d = std::numeric_limits<float>::max();
                    for (int i=0; i<3; ++i) {
                        float d = color_dist(p, c[i], true);
                        if (d < min_d) { min_d = d; best_idx = i; }
                    }
                }
                int bit_offset = (col < 4 ? 0 : 0x20) + 2 * (r * 4 + (col & 3));
                set_bit(data, bit_offset, best_idx & 1);
                set_bit(data, bit_offset+1, (best_idx >> 1) & 1);
            }
        } else {
            for (int block_idx=0; block_idx<2; ++block_idx) {
                uint32_t palette[4][4];
                if (block_idx==0) for(int j=0;j<4;++j) for(int k=0;k<4;++k) palette[j][k] = (c[0].components[k]*(3-j) + c[1].components[k]*j+1)/3;
                else for(int j=0;j<4;++j) for(int k=0;k<4;++k) palette[j][k] = (c[2].components[k]*(3-j) + c[1].components[k]*j+1)/3;
                
                int bit_stream_offset = dword_4B3B90[block_idx];
                for (int r=0; r<4; ++r) for (int col=0; col<4; ++col) {
                    Color p = block_pixels[r][block_idx*4+col];
                    int best_idx = 0;
                    float min_d = std::numeric_limits<float>::max();
                    for (int i=0; i<4; ++i) {
                        Color pal_c = {(uint8_t)palette[i][0], (uint8_t)palette[i][1], (uint8_t)palette[i][2], (uint8_t)palette[i][3]};
                        float d = color_dist(p, pal_c, true);
                        if (d < min_d) { min_d = d; best_idx = i; }
                    }
                    int bit_offset = r * 8 + col * 2;
                    set_bit(data, bit_stream_offset + bit_offset, best_idx & 1);
                    set_bit(data, bit_stream_offset + bit_offset+1, (best_idx >> 1) & 1);
                }
            }
        }
        data[3] = a8;
        float err = compute_block_error(block_pixels, data);
        if (err < best_error) { best_error = err; best_data = data; }
    }
    
    return best_data;
}

// Compress image
void compress_image(const unsigned char* image_data, int width, int height, std::vector<char>& compressed) {
    int blocks_x = (width + 7) / 8;
    int blocks_y = (height + 3) / 4;
    for (int by = 0; by < blocks_y; ++by) {
        for (int bx = 0; bx < blocks_x; ++bx) {
            Color block[4][8];
            extract_block(image_data, width, height, bx, by, block);
            std::vector<uint32_t> block_data = compress_block(block);
            const char* ptr = reinterpret_cast<const char*>(block_data.data());
            compressed.insert(compressed.end(), ptr, ptr + 16);
        }
    }
}