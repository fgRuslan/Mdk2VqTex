#pragma once

#include <string>

void dump_flags(char* texture, std::string output_file);
void load_flags(char* texture, std::string input_file);
bool is_compression_disabled_in_ini(std::string ini_path);