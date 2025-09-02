#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#pragma pack(push,1)
struct TextureHeader {
    uint32_t magic;
    uint32_t version;
    uint32_t width;
    uint32_t height;
    uint32_t flags;
    uint32_t data_size;
};
#pragma pack(pop)

// --- Lookup Tables (Values from original code) ---
const uint32_t dword_4B3B80[] = { 0x40, 0x5E };
const uint32_t dword_4B3B88[] = { 0, 0x20 };
const uint32_t dword_4B3B90[] = { 0, 0x20 };
// ===============================================================================================
// UTILITY FUNCTIONS
// ===============================================================================================

// Helper to safely read from the compressed data buffer
inline uint32_t get_bit(const std::vector<uint32_t>& data, int bit_index) {
    if ((bit_index >> 5) >= data.size()) return 0;
    return (data[bit_index >> 5] >> (bit_index & 0x1F)) & 1;
}

// ===============================================================================================
// PORTED DECOMPRESSION FUNCTIONS
// ===============================================================================================

// Forward declarations
void sub_45B6D0(int width, int x, int y, uint32_t* image_buffer, const std::vector<uint32_t>& pixel_data, uint32_t a8);
void sub_45B930(int width, int x, int y, uint32_t* image_buffer, const std::vector<uint32_t>& pixel_data);
void sub_45BB00(int width, int x, int y, uint32_t* image_buffer, const std::vector<uint32_t>& pixel_data);
void sub_45BCE0(int width, int x, int y, uint32_t* image_buffer, const std::vector<uint32_t>& pixel_data, uint32_t a8);


// Corresponds to sub_45B530
// Generates a palette based on the exact reverse-order interpolation logic from the original game.
void generate_palette_from_ida(uint32_t* palette, const uint32_t* c1, const uint32_t* c2, bool is_3_color_mode) {
    if (!is_3_color_mode) {
        // --- 4-Color Opaque Mode ---
        // CRITICAL: The interpolation is from c2 to c1, not the other way around.
        // Palette[0] = c2
        // Palette[1] = 1/3 c1 + 2/3 c2
        // Palette[2] = 2/3 c1 + 1/3 c2
        // Palette[3] = c1
        for (int i = 0; i < 3; ++i) { // RGB components
            palette[0 * 4 + i] = c2[i];
            palette[1 * 4 + i] = (1 * c1[i] + 2 * c2[i] + 1) / 3;
            palette[2 * 4 + i] = (2 * c1[i] + 1 * c2[i] + 1) / 3;
            palette[3 * 4 + i] = c1[i];
        }
        for (int i = 0; i < 4; ++i) {
            palette[i * 4 + 3] = 255; // Set full alpha for all colors
        }
    } else {
        // --- 3-Color + 1-Bit Alpha Mode ---
        // CRITICAL: The interpolation is also from c2 to c1 here.
        // Palette[0] = c2
        // Palette[1] = 1/2 c1 + 1/2 c2
        // Palette[2] = c1
        // Palette[3] = Transparent
        for (int i = 0; i < 3; ++i) { // RGB components
            palette[0 * 4 + i] = c2[i];
            palette[1 * 4 + i] = (c1[i] + c2[i]) / 2; // Note: Original uses integer division, no +1 rounding
            palette[2 * 4 + i] = c1[i];
        }
        palette[0 * 4 + 3] = 255;
        palette[1 * 4 + 3] = 255;
        palette[2 * 4 + 3] = 255;

        // Set the fourth color to transparent black
        palette[3 * 4 + 0] = 0;
        palette[3 * 4 + 1] = 0;
        palette[3 * 4 + 2] = 0;
        palette[3 * 4 + 3] = 0;
    }
}

// Corresponds to sub_45B630
// Generates an 8-color palette by interpolating between two RGB colors.
void generate_palette_mode1(uint32_t* palette, const uint32_t* colors1, const uint32_t* colors2) {
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 3; ++j) { // Original only processes 3 channels (RGB)
            palette[i * 4 + j] = (colors1[j] * (7 - i) + colors2[j] * i + 3) / 7;
        }
        palette[i * 4 + 3] = 0xFF; // Set alpha to full opacity
    }
}


// Corresponds to sub_45B6D0
void sub_45B6D0(int width, int x, int y, uint32_t* image_buffer, const std::vector<uint32_t>& pixel_data, uint32_t a8) {
    bool v27 = (a8 >> 28) & 1; // This is the 'a1' flag for the palette function

    for (int block_idx = 0; block_idx < 2; ++block_idx) {
        int block_offset_x = block_idx * 4;
        uint32_t c1[4] = {0}, c2[4] = {0}; // c1 is a4, c2 is a3

        uint32_t color_base_offset = dword_4B3B80[block_idx];
        bool v11 = (block_idx == 0) ? ((a8 >> 29) & 1) : ((a8 >> 30) & 1);

        for (int i = 0; i < 4; ++i) { // RGBA channels
            uint32_t c1_comp_offset = color_base_offset;
            uint32_t c2_comp_offset = color_base_offset + 15;
            for (int bit_idx = 0; bit_idx < 5; ++bit_idx) {
                c1[i] |= get_bit(pixel_data, c1_comp_offset + i * 5 + bit_idx) << (3 + bit_idx);
                c2[i] |= get_bit(pixel_data, c2_comp_offset + i * 5 + bit_idx) << (3 + bit_idx);
            }
            c1[i] |= (c1[i] >> 5);
            c2[i] |= (c2[i] >> 5);
        }

        c1[3] = (c1[3] & 0xF8) | (4 * v11);
        if (!v27) {
             c2[3] = (c2[3] & 0xF8) | (4 * (v11 ^ get_bit(pixel_data, dword_4B3B88[block_idx] + 1)));
        }
        
        uint32_t palette[4 * 4];
        
        generate_palette_from_ida(palette, c1, c2, v27);

        uint32_t* block_ptr = image_buffer + (y * width) + x + block_offset_x;
        int bit_stream_offset = dword_4B3B88[block_idx];

        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col) {
                int bit_offset = row * 8 + col * 2;
                uint32_t index_bit0 = get_bit(pixel_data, bit_stream_offset + bit_offset);
                uint32_t index_bit1 = get_bit(pixel_data, bit_stream_offset + bit_offset + 1);
                uint32_t palette_index = index_bit0 | (index_bit1 << 1);

                uint32_t r = palette[palette_index * 4 + 0];
                uint32_t g = palette[palette_index * 4 + 1];
                uint32_t b = palette[palette_index * 4 + 2];
                uint32_t a = palette[palette_index * 4 + 3];

                if ((y + row < 1024) && (x + block_offset_x + col < 1024)) {
                   block_ptr[row * width + col] = (a << 24) | (b << 16) | (g << 8) | r;
                }
            }
        }
    }
}

// Corresponds to sub_45B930
void sub_45B930(int width, int x, int y, uint32_t* image_buffer, const std::vector<uint32_t>& pixel_data) {
    uint32_t c1[3] = {0}, c2[3] = {0};
    
    for (int i = 0; i < 3; ++i) { // RGB
        int component_offset = i * 5;
        for (int bit_idx = 0; bit_idx < 5; ++bit_idx) {
            c1[i] |= get_bit(pixel_data, 96 + component_offset + bit_idx) << (3 + bit_idx);
            c2[i] |= get_bit(pixel_data, 111 + component_offset + bit_idx) << (3 + bit_idx);
        }
        c1[i] |= c1[i] >> 5;
        c2[i] |= c2[i] >> 5;
    }

    uint32_t palette[8 * 4]; // 8 colors, 4 channels
    generate_palette_mode1(palette, c1, c2);
    
    uint32_t* block_ptr = image_buffer + y * width + x;
    
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 8; ++col) {
            // Correct bit offset calculation based on column position
            int bit_offset = (col < 4 ? 0 : 0x30) + 3 * ((row * 4) + (col & 3));
            
            uint32_t idx0 = get_bit(pixel_data, bit_offset);
            uint32_t idx1 = get_bit(pixel_data, bit_offset + 1);
            uint32_t idx2 = get_bit(pixel_data, bit_offset + 2);
            uint32_t palette_index = idx0 | (idx1 << 1) | (idx2 << 2);
            
            uint32_t r = palette[palette_index * 4 + 0];
            uint32_t g = palette[palette_index * 4 + 1];
            uint32_t b = palette[palette_index * 4 + 2];
            uint32_t a = palette[palette_index * 4 + 3];
            
            if (y + row < 1024 && x + col < 1024) {
                block_ptr[row * width + col] = (a << 24) | (b << 16) | (g << 8) | r;
            }
        }
    }
}

// Corresponds to sub_45BB00
void sub_45BB00(int width, int x, int y, uint32_t* image_buffer, const std::vector<uint32_t>& pixel_data) {
    uint32_t colors[4][4] = {{0}};

    int base_bit_index = 79;
    // Outer loop iterates through components (R, G, B)
    for (int comp_idx = 0; comp_idx < 3; ++comp_idx) {
        // Inner loop builds the 5 bits for the current component for each color.
        for (int bit_pos = 0; bit_pos < 5; ++bit_pos) {
            int current_bit_index = base_bit_index + bit_pos;
            int dest_bit = 3 + bit_pos;

            // Read the bit for the current component (comp_idx) for each of the 4 colors.
            // The source bit locations are interleaved in the data stream.
            colors[0][comp_idx] |= get_bit(pixel_data, current_bit_index - 15) << dest_bit;
            colors[1][comp_idx] |= get_bit(pixel_data, current_bit_index)      << dest_bit;
            colors[2][comp_idx] |= get_bit(pixel_data, current_bit_index + 15) << dest_bit;
            colors[3][comp_idx] |= get_bit(pixel_data, current_bit_index + 30) << dest_bit;
        }
        base_bit_index += 5; // Move base index to the start of the next component block.
    }

    // Finalize the colors: apply the ">> 5" replication and set alpha.
    for (int c_idx = 0; c_idx < 4; ++c_idx) {
        for (int comp_idx = 0; comp_idx < 3; ++comp_idx) {
            colors[c_idx][comp_idx] |= colors[c_idx][comp_idx] >> 5;
        }
        colors[c_idx][3] = 255; // Alpha is always opaque.
    }

    uint32_t* block_ptr = image_buffer + y * width + x;

    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 8; ++col) {
            int bit_offset = (col < 4 ? 0 : 0x20) + 2 * ((row * 4) + (col & 3));
            uint32_t idx0 = get_bit(pixel_data, bit_offset);
            uint32_t idx1 = get_bit(pixel_data, bit_offset + 1);
            uint32_t palette_index = idx0 | (idx1 << 1);

            uint32_t r = colors[palette_index][0];
            uint32_t g = colors[palette_index][1];
            uint32_t b = colors[palette_index][2];
            uint32_t a = colors[palette_index][3];

            if (y + row < 1024 && x + col < 1024) {
                block_ptr[row * width + col] = (a << 24) | (b << 16) | (g << 8) | r;
            }
        }
    }
}

// Corresponds to sub_45BCE0
void sub_45BCE0(int width, int x, int y, uint32_t* image_buffer, const std::vector<uint32_t>& pixel_data, uint32_t a8) {
    bool v58 = (a8 >> 28) & 1;
    uint32_t c[3][4] = {{0}};

    // This function also reads color data component-first, like sub_45BB00.
    // It reads the R for all 3 colors, then G for all 3, then B for all 3.
    int base_bit_index = 79;
    // Outer loop iterates through components (R, G, B)
    for (int comp_idx = 0; comp_idx < 3; ++comp_idx) {
        // Inner loop builds the 5 bits for the current component for each of the 3 colors.
        for (int bit_pos = 0; bit_pos < 5; ++bit_pos) {
            int current_bit_index = base_bit_index + bit_pos;
            int dest_bit = 3 + bit_pos;

            // Read the bit for the current component (comp_idx) for each of the 3 colors.
            c[0][comp_idx] |= get_bit(pixel_data, current_bit_index - 15) << dest_bit;
            c[1][comp_idx] |= get_bit(pixel_data, current_bit_index)      << dest_bit;
            c[2][comp_idx] |= get_bit(pixel_data, current_bit_index + 15) << dest_bit;
        }
        base_bit_index += 5; // Move to the start of the next component's data.
    }

    // Finalize the colors: apply the ">> 5" replication and set alpha.
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            c[i][j] |= c[i][j] >> 5;
        }
        c[i][3] = 255; // Set alpha to full opacity.
    }

    // The rest of the logic for both the interpolated (v58) and non-interpolated
    // paths was correct and can remain unchanged.
    if (v58) {
        for (int block_idx = 0; block_idx < 2; ++block_idx) {
            uint32_t palette[4][4];
            uint32_t* cA = (block_idx == 0) ? c[0] : c[1];
            uint32_t* cB = (block_idx == 0) ? c[1] : c[2];

            for(int j=0; j<4; ++j) {
                for(int k=0; k<4; ++k) {
                    palette[j][k] = (cA[k] * (3-j) + cB[k] * j + 1) / 3;
                }
            }
            
            uint32_t* block_ptr = image_buffer + (y * width) + x + (block_idx * 4);
            int bit_stream_offset = dword_4B3B90[block_idx];

            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                     int bit_offset = row*8 + col*2;
                     uint32_t idx0 = get_bit(pixel_data, bit_stream_offset + bit_offset);
                     uint32_t idx1 = get_bit(pixel_data, bit_stream_offset + bit_offset + 1);
                     uint32_t pal_idx = idx0 | (idx1 << 1);

                     uint32_t r = palette[pal_idx][0];
                     uint32_t g = palette[pal_idx][1];
                     uint32_t b = palette[pal_idx][2];
                     uint32_t a = palette[pal_idx][3];
                     if (y + row < 1024 && x + col + block_idx*4 < 1024) {
                         block_ptr[row*width+col] = (a << 24) | (b << 16) | (g << 8) | r;
                     }
                }
            }
        }
    } else {
        uint32_t* block_ptr = image_buffer + y * width + x;
        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 8; ++col) {
                int bit_offset = (col < 4 ? 0 : 0x20) + 2 * ((row * 4) + (col & 3));
                uint32_t idx0 = get_bit(pixel_data, bit_offset);
                uint32_t idx1 = get_bit(pixel_data, bit_offset + 1);
                uint32_t pal_idx = idx0 | (idx1 << 1);
                
                uint32_t r = (pal_idx < 3) ? c[pal_idx][0] : 0;
                uint32_t g = (pal_idx < 3) ? c[pal_idx][1] : 0;
                uint32_t b = (pal_idx < 3) ? c[pal_idx][2] : 0;
                uint32_t a = (pal_idx < 3) ? c[pal_idx][3] : 255;

                if (y + row < 1024 && x + col < 1024) {
                    block_ptr[row * width + col] = (a << 24) | (b << 16) | (g << 8) | r;
                }
            }
        }
    }
}


// Corresponds to sub_45C130, the main dispatcher
void decode_block(int width, int x, int y, uint32_t* image_buffer, const std::vector<uint32_t>& pixel_data) {
    uint32_t a8 = pixel_data[3];

    if ((a8 & 0x80000000) != 0) {
        sub_45B6D0(width, x, y, image_buffer, pixel_data, a8);
    } else if (((a8 >> 30) & 1) == 0) {
        sub_45B930(width, x, y, image_buffer, pixel_data);
    } else if (((a8 >> 29) & 1) != 0) {
        sub_45BCE0(width, x, y, image_buffer, pixel_data, a8);
    } else {
        sub_45BB00(width, x, y, image_buffer, pixel_data);
    }
}

// Corresponds to ProcessTextureBlocks
void decompress_image(uint32_t* decompressed_output, const char* compressed_input, int width, int height) {
    const uint32_t* data_ptr = reinterpret_cast<const uint32_t*>(compressed_input);

    for (int y = 0; y < height; y += 4) {
        for (int x = 0; x < width; x += 8) {
            std::vector<uint32_t> block_data(4);
            // The original code passes the whole block as arguments, so we load it here.
            block_data[0] = data_ptr[0]; // Corresponds to a5
            block_data[1] = data_ptr[1]; // Corresponds to a6
            block_data[2] = data_ptr[2]; // Corresponds to a7
            block_data[3] = data_ptr[3]; // Corresponds to a8
            data_ptr += 4;
            decode_block(width, x, y, decompressed_output, block_data);
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_texture_file> <output_png_file>" << std::endl;
        return 1;
    }
    std::string input_path = argv[1];
    std::string output_path = argv[2];
    std::ifstream file(input_path, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Cannot open input file: " << input_path << std::endl;
        return 1;
    }
    file.seekg(0, std::ios::end);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);
    std::vector<char> buffer(size);
    if (!file.read(buffer.data(), size)) {
        std::cerr << "Error: Failed to read from input file." << std::endl;
        return 1;
    }
    if (buffer.size() < sizeof(TextureHeader)) {
        std::cerr << "Error: File is too small to contain a valid header." << std::endl;
        return 1;
    }
    const TextureHeader* header = reinterpret_cast<const TextureHeader*>(buffer.data());
    int width = header->width;
    int height = header->height;
    if (width <= 0 || width > 4096 || height <= 0 || height > 4096) {
        std::cerr << "Error: Invalid texture dimensions in header (" << width << "x" << height << ")." << std::endl;
        std::cerr << "You may need to adjust the TextureHeader struct or file offset." << std::endl;
        return 1;
    }
    std::cout << "Texture Info:" << std::endl;
    std::cout << "  Dimensions: " << width << "x" << height << std::endl;
    std::vector<uint32_t> decompressed_image_buffer(width * height);
    const char* pixel_data_ptr = buffer.data() + sizeof(TextureHeader);
    decompress_image(decompressed_image_buffer.data(), pixel_data_ptr, width, height);
    for (uint32_t& pixel : decompressed_image_buffer) {
        uint32_t a = (pixel >> 24) & 0xFF;
        uint32_t r = (pixel >> 0)  & 0xFF;
        uint32_t g = (pixel >> 8)  & 0xFF;
        uint32_t b = (pixel >> 16) & 0xFF;
        pixel = (a << 24) | (b << 16) | (g << 8) | r;
    }
    int channels = 4;
    if (stbi_write_png(output_path.c_str(), width, height, channels, decompressed_image_buffer.data(), width * channels)) {
        std::cout << "Successfully decompressed and saved texture to " << output_path << std::endl;
    } else {
        std::cerr << "Error: Failed to write output PNG file." << std::endl;
        return 1;
    }

    return 0;
}