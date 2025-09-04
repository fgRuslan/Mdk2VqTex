#pragma once

#include <cstdint>

// ===============================================================================================
// UTILITY FUNCTIONS
// ===============================================================================================

// Helper to safely read from the compressed data buffer
inline uint32_t get_bit(const std::vector<uint32_t>& data, int bit_index) {
    if ((bit_index >> 5) >= data.size()) return 0;
    return (data[bit_index >> 5] >> (bit_index & 0x1F)) & 1;
}

void generate_palette_from_ida(uint32_t* palette, const uint32_t* c1, const uint32_t* c2, bool is_3_color_mode);
void generate_palette_mode1(uint32_t* palette, const uint32_t* colors1, const uint32_t* colors2);

void decode_block(int width, int x, int y, uint32_t* image_buffer, const std::vector<uint32_t>& pixel_data);
void decompress_image(uint32_t* decompressed_output, const char* compressed_input, int width, int height);