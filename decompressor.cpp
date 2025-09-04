#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cctype>

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
    }
    else {
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
void generate_palette_mode1(uint32_t* palette, const uint32_t* colors1, const uint32_t* colors2)
{
    // This function generates 7 interpolated colors and one transparent color.
    // The loop iterates 7 times to generate the main part of the palette.
    for (int i = 0; i < 7; ++i) {
        uint32_t row = i;         // Corresponds to the weight of colors1
        uint32_t weight = 6 - i;  // Corresponds to the weight of colors2

        uint32_t* output_ptr = palette + (i * 4); // Point to the start of the current color in the palette

        // Interpolate the R, G, B components
        for (int comp = 0; comp < 3; ++comp) {
            uint32_t term1 = colors1[comp] * weight;
            uint32_t term2 = colors2[comp] * row;
            output_ptr[comp] = (term1 + term2 + 2) / 6;
        }

        output_ptr[3] = 0xFF;  // Set alpha channel to opaque
    }

    // Set the 8th and final color to transparent black
    palette[28] = 0;
    palette[29] = 0;
    palette[30] = 0;
    palette[31] = 0;
}


// Corresponds to sub_45B6D0
void sub_45B6D0(int width, int x, int y, uint32_t* image_buffer, const std::vector<uint32_t>& pixel_data, uint32_t a8) {
    bool v27 = (a8 >> 28) & 1; // This is the 'a1' flag for the palette function

    for (int block_idx = 0; block_idx < 2; ++block_idx) {
        int block_offset_x = block_idx * 4;
        uint32_t c1[4] = { 0 }, c2[4] = { 0 }; // c1 is a4, c2 is a3

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
    uint32_t c1[3] = { 0 }, c2[3] = { 0 };

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
    uint32_t colors[4][4] = { {0} };

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
            colors[1][comp_idx] |= get_bit(pixel_data, current_bit_index) << dest_bit;
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
    uint32_t c[3][4] = { {0} };

    // This function reads color data component-first.
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
            c[1][comp_idx] |= get_bit(pixel_data, current_bit_index) << dest_bit;
            c[2][comp_idx] |= get_bit(pixel_data, current_bit_index + 15) << dest_bit;
        }
        base_bit_index += 5; // Move to the start of the next component's data.
    }

    // Read the 5-bit Alpha for all 3 colors, matching the decompiled version's logic.
    for (int bit_pos = 0; bit_pos < 5; ++bit_pos) {
        int dest_bit = 3 + bit_pos;
        int base_idx = 13 + bit_pos;
        c[0][3] |= get_bit(pixel_data, base_idx + 96) << dest_bit;
        c[1][3] |= get_bit(pixel_data, base_idx + 101) << dest_bit;
        c[2][3] |= get_bit(pixel_data, base_idx + 106) << dest_bit;
    }

    // Finalize the colors: apply the ">> 5" replication to all RGBA components.
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) { // Loop to 4 to include Alpha
            c[i][j] |= c[i][j] >> 5;
        }
    }

    // Interpolated path
    if (v58) {
        for (int block_idx = 0; block_idx < 2; ++block_idx) {
            uint32_t palette[4][4];

            if (block_idx == 0) {
                // First block interpolates from c[0] to c[1]
                for (int j = 0; j < 4; ++j) {
                    for (int k = 0; k < 4; ++k) {
                        palette[j][k] = (c[0][k] * (3 - j) + c[1][k] * j + 1) / 3;
                    }
                }
            }
            else {
                // Second block interpolates from c[2] to c[1]
                for (int j = 0; j < 4; ++j) {
                    for (int k = 0; k < 4; ++k) {
                        palette[j][k] = (c[2][k] * (3 - j) + c[1][k] * j + 1) / 3;
                    }
                }
            }

            uint32_t* block_ptr = image_buffer + (y * width) + x + (block_idx * 4);
            int bit_stream_offset = dword_4B3B90[block_idx];

            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    int bit_offset = row * 8 + col * 2;
                    uint32_t idx0 = get_bit(pixel_data, bit_stream_offset + bit_offset);
                    uint32_t idx1 = get_bit(pixel_data, bit_stream_offset + bit_offset + 1);
                    uint32_t pal_idx = idx0 | (idx1 << 1);

                    uint32_t r = palette[pal_idx][0];
                    uint32_t g = palette[pal_idx][1];
                    uint32_t b = palette[pal_idx][2];
                    uint32_t a = palette[pal_idx][3];
                    if (y + row < 1024 && x + col + block_idx * 4 < 1024) {
                        block_ptr[row * width + col] = (a << 24) | (b << 16) | (g << 8) | r;
                    }
                }
            }
        }
    }
    else {
        // Non-interpolated path
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
                uint32_t a = (pal_idx < 3) ? c[pal_idx][3] : 0;

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
    }
    else if (((a8 >> 30) & 1) == 0) {
        sub_45B930(width, x, y, image_buffer, pixel_data);
    }
    else if (((a8 >> 29) & 1) != 0) {
        sub_45BCE0(width, x, y, image_buffer, pixel_data, a8);
    }
    else {
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