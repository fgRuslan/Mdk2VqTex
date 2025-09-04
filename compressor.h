#pragma once
#include <vector>

void init_quant_tables();
void compress_image(const unsigned char* image_data, int width, int height, std::vector<char>& compressed);