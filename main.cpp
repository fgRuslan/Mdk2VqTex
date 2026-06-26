#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cctype>
#include <string>

#include "decompressor.h"
#include "compressor.h"
#include "flag_processor.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "stb_image_resize.h"

//#include "SimpleIni.h"

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

bool endsWith(std::string const& fullString, std::string const& ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    }
    else {
        return false;
    }
}

int decompress(std::string input_path, std::string output_path, bool do_flip, bool decompress_all_mips)
{
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
    if (width <= 0 || height <= 0) {
        std::cerr << "Error: Invalid texture dimensions in header (" << width << "x" << height << ")." << std::endl;
    }
    std::cout << "Texture Info:" << std::endl;
    std::cout << "  Dimensions: " << width << "x" << height << std::endl;

    const uint32_t* metadata_u32 = reinterpret_cast<const uint32_t*>(buffer.data() + sizeof(TextureHeader));
    bool is_compressed = (metadata_u32[3] & 0x20) != 0;

    if (!is_compressed) {
        std::cout << "The texture is UNCOMPRESSED. Reading raw pixel data..." << std::endl;

        int channels = header->flags;

        int data_offset = sizeof(TextureHeader) + 9 * sizeof(uint32_t);

        std::cout << data_offset + width * height * channels << std::endl;
        std::cout << buffer.size() << std::endl;

        std::vector<uint8_t> rgba_pixels(width * height * 4);
        const uint8_t* raw_data = reinterpret_cast<const uint8_t*>(buffer.data() + data_offset);

        for (int i = 0; i < width * height; ++i) {
            rgba_pixels[i * 4 + 0] = raw_data[i * channels + 2];
            rgba_pixels[i * 4 + 1] = raw_data[i * channels + 1];
            rgba_pixels[i * 4 + 2] = raw_data[i * channels + 0];
            rgba_pixels[i * 4 + 3] = (channels == 4) ? raw_data[i * channels + 3] : 255;
        }

        if (do_flip) {
            std::cout << "Performing a flip on Y axis..." << std::endl;
            for (int y = 0; y < height / 2; ++y) {
                uint8_t* row1 = rgba_pixels.data() + y * width * 4;
                uint8_t* row2 = rgba_pixels.data() + (height - 1 - y) * width * 4;
                for (int x = 0; x < width * 4; ++x) {
                    std::swap(row1[x], row2[x]);
                }
            }
        }

        if (stbi_write_png(output_path.c_str(), width, height, 4, rgba_pixels.data(), width * 4)) {
            std::cout << "Successfully saved uncompressed texture to " << output_path << std::endl;
        }
        else {
            std::cerr << "Error: Failed to write output PNG file." << std::endl;
            return 1;
        }

        std::string base_name = output_path.substr(0, output_path.find_last_of("."));
        dump_flags((char*)buffer.data(), base_name + ".ini");

        return 0;
    }

    if (decompress_all_mips) {
        std::cout << "Decompressing all mipmaps..." << std::endl;

        uint32_t metadata_u32_count = 9;//9 base fields
        //scan starting from index 9 (the last base field) to find valid offsets
        for (int i = 9; i < (size / sizeof(uint32_t)); ++i) {
            uint32_t offset = metadata_u32[i];
            if (offset > 0 && offset < size) {
                metadata_u32_count = i + 1;
            } else {
                break;
            }
        }

        // The compressor writes offsets in reverse order (smallest mip first).
        // We read them into a temporary vector.
        //TODO: check if it's the right order
        std::vector<uint32_t> mip_offsets_rev;
        for (int i = 9; i < metadata_u32_count; ++i) {
            uint32_t offset = metadata_u32[i];
            if (offset > 0 && offset < size) {
                mip_offsets_rev.push_back(offset);
            }
            else {
                break; // Stop at the first invalid or zero offset
            }
        }

        if (mip_offsets_rev.empty()) {
            std::cerr << "Error: No mipmap offset data found in the header. Cannot extract all mips." << std::endl;
            return 1;
        }

        std::vector<uint32_t> mip_offsets = mip_offsets_rev;
        std::reverse(mip_offsets.begin(), mip_offsets.end());

        std::string base_name;
        std::string extension;
        size_t dot_pos = output_path.find_last_of(".");
        if (dot_pos != std::string::npos) {
            base_name = output_path.substr(0, dot_pos);
            extension = output_path.substr(dot_pos);
        }
        else {
            base_name = output_path;
            extension = ".png";
        }

        for (size_t i = 0; i < mip_offsets.size(); ++i) {
            int mip_w = std::max(1, (int)header->width >> i);
            int mip_h = std::max(1, (int)header->height >> i);

            std::cout << "  Decompressing Mip " << i << " (" << mip_w << "x" << mip_h << ")..." << std::endl;

            std::vector<uint32_t> decompressed_image_buffer(mip_w * mip_h);
            const char* pixel_data_ptr = buffer.data() + mip_offsets[i];
            decompress_image(decompressed_image_buffer.data(), pixel_data_ptr, mip_w, mip_h);

            if (do_flip) {
                for (int y = 0; y < mip_h / 2; ++y) {
                    auto* row1 = decompressed_image_buffer.data() + y * mip_w;
                    auto* row2 = decompressed_image_buffer.data() + (mip_h - 1 - y) * mip_w;
                    std::swap_ranges(row1, row1 + mip_w, row2);
                }
            }

            for (uint32_t& pixel : decompressed_image_buffer) {
                uint32_t a = (pixel >> 24) & 0xFF;
                uint32_t r = (pixel >> 16) & 0xFF;
                uint32_t g = (pixel >> 8) & 0xFF;
                uint32_t b = (pixel >> 0) & 0xFF;
                pixel = (a << 24) | (b << 16) | (g << 8) | r;
            }

            std::string mip_output_path = base_name + "_mip" + std::to_string(i) + extension;
            int channels = 4;
            if (stbi_write_png(mip_output_path.c_str(), mip_w, mip_h, channels, decompressed_image_buffer.data(), mip_w * channels)) {
                std::cout << "    Successfully saved to " << mip_output_path << std::endl;
            }
            else {
                std::cerr << "    Error: Failed to write output PNG file: " << mip_output_path << std::endl;
            }
        }

    }
    else {
        //here we decompress only the largest mip
		//TODO: probably I shouldn't repeat myself and just use the same code as in the loop above, but stop after mip0
        std::vector<uint32_t> decompressed_image_buffer(width * height);

        const uint32_t* metadata_u32 = reinterpret_cast<const uint32_t*>(buffer.data() + sizeof(TextureHeader));

        uint32_t metadata_u32_count = 9;//9 base fields
        //scan starting from index 9 (the last base field) to find valid offsets
        for (int i = 9; i < (size / sizeof(uint32_t)); ++i) {
            uint32_t offset = metadata_u32[i];
            if (offset > 0 && offset < size) {
                metadata_u32_count = i + 1;
            }
            else {
                break;
            }
        }

        std::vector<uint32_t> mip_offsets_rev;
        for (int i = 9; i < metadata_u32_count; ++i) {
            uint32_t offset = metadata_u32[i];
            if (offset > 0 && offset < size) {
                mip_offsets_rev.push_back(offset);
            }
            else {
                break; // Stop at the first invalid or zero offset
            }
        }

        if (mip_offsets_rev.empty()) {
            std::cerr << "Error: No mipmap offset data found in the header. Cannot extract all mips." << std::endl;
            return 1;
        }

        std::vector<uint32_t> mip_offsets = mip_offsets_rev;
        std::reverse(mip_offsets.begin(), mip_offsets.end());

        const char* pixel_data_ptr = buffer.data() + mip_offsets[0];
        decompress_image(decompressed_image_buffer.data(), pixel_data_ptr, width, height);

        if (do_flip)
        {
            std::cout << "Performing a flip on Y axis..." << std::endl;
            for (int y = 0; y < height / 2; ++y)
            {
                auto* row1 = decompressed_image_buffer.data() + y * width;
                auto* row2 = decompressed_image_buffer.data() + (height - 1 - y) * width;
                std::swap_ranges(row1, row1 + width, row2);
            }
        }

        for (uint32_t& pixel : decompressed_image_buffer) {
            uint32_t a = (pixel >> 24) & 0xFF;
            uint32_t r = (pixel >> 16) & 0xFF;
            uint32_t g = (pixel >> 8) & 0xFF;
            uint32_t b = (pixel >> 0) & 0xFF;
            pixel = (a << 24) | (b << 16) | (g << 8) | r;
        }

        int channels = 4;
        if (stbi_write_png(output_path.c_str(), width, height, channels, decompressed_image_buffer.data(), width * channels)) {
            std::cout << "Successfully decompressed and saved texture to " << output_path << std::endl;
        }
        else {
            std::cerr << "Error: Failed to write output PNG file." << std::endl;
            return 1;
        }
    }

    std::string base_name;
    size_t dot_pos = output_path.find_last_of(".");
    if (dot_pos != std::string::npos) {
        base_name = output_path.substr(0, dot_pos);
    }
    else {
        base_name = output_path;
    }

    dump_flags(buffer.data(), base_name + ".ini");

    return 0;
}

int compress(std::string input_path, std::string output_path, bool do_flip, bool uncompressed)
{
    int w, h, channels;
	uint32_t contains_alpha = 0x0;
	
    unsigned char* image_data = stbi_load(input_path.c_str(), &w, &h, &channels, 4);
    if (!image_data) {
        std::cerr << "Error: Cannot open or read input PNG: " << input_path << std::endl;
        return 1;
    }
    std::cout << "Texture Info:" << std::endl;
    std::cout << "  Dimensions: " << w << "x" << h << std::endl;

    if (do_flip) {
        std::cout << "Performing a flip on Y axis..." << std::endl;
        for (int y = 0; y < h / 2; ++y) {
            unsigned char* row1 = image_data + y * w * 4;
            unsigned char* row2 = image_data + (h - 1 - y) * w * 4;
            for (int x = 0; x < w * 4; ++x) {
                std::swap(row1[x], row2[x]);
            }
        }
    }

    std::string base_name = input_path.substr(0, input_path.find_last_of("."));
    std::string ini_path = base_name + ".ini";

    bool is_uncompressed = is_compression_disabled_in_ini(ini_path);

    if (is_uncompressed) {
        std::cout << "Saving in UNCOMPRESSED mode..." << std::endl;

        int tex_channels = (channels == 3) ? 3 : 4;
        uint32_t contains_alpha = 0x0;

        if (tex_channels == 4) {
            for (int i = 0; i < w * h * 4; i += 4) {
                if (image_data[i + 3] < 255) {
                    contains_alpha = 0x01;
                    break;
                }
            }
        }

        TextureHeader header;
        header.width = w;
        header.height = h;
        header.flags = tex_channels;
        header.data_size = w * h * tex_channels;

        uint32_t metadata_size = 9 * sizeof(uint32_t);
        std::vector<char> metadata(metadata_size, 0);
        uint32_t* metadata_u32 = reinterpret_cast<uint32_t*>(metadata.data());

        metadata_u32[2] = contains_alpha;
        metadata_u32[3] = 0x00;//no compression
        metadata_u32[5] = 0x54455843;
        metadata_u32[7] = w;
        metadata_u32[8] = h;

        std::vector<char> file_buffer;
        file_buffer.reserve(sizeof(TextureHeader) + metadata_size + header.data_size);

        const char* header_ptr = reinterpret_cast<const char*>(&header);
        file_buffer.insert(file_buffer.end(), header_ptr, header_ptr + sizeof(header));
        file_buffer.insert(file_buffer.end(), metadata.data(), metadata.data() + metadata.size());

        std::vector<char> raw_pixels(w * h * tex_channels);
        for (int i = 0; i < w * h; ++i) {
            raw_pixels[i * tex_channels + 0] = image_data[i * 4 + 2];//b
            raw_pixels[i * tex_channels + 1] = image_data[i * 4 + 1];//g
            raw_pixels[i * tex_channels + 2] = image_data[i * 4 + 0];//r
            if (tex_channels == 4) {
                raw_pixels[i * tex_channels + 3] = image_data[i * 4 + 3];//a
            }
        }
        file_buffer.insert(file_buffer.end(), raw_pixels.begin(), raw_pixels.end());

        load_flags(file_buffer.data(), ini_path);

        std::ofstream out_file(output_path, std::ios::binary);
        out_file.write(file_buffer.data(), file_buffer.size());

        std::cout << "Successfully saved uncompressed texture to " << output_path << std::endl;
        stbi_image_free(image_data);
        return 0;
    }

    init_quant_tables();

    // --- Start of Mipmap Generation and Compression ---

    std::vector<std::vector<char>> mip_data_levels;
    std::vector<uint32_t> mip_offsets;
    std::vector<uint32_t> mip_sizes;
    uint32_t total_compressed_size = 0;

    int current_w = w;
    int current_h = h;
    unsigned char* current_image_data = image_data;
    unsigned char* prev_image_data = nullptr;

    std::cout << "Generating and compressing mipmaps..." << std::endl;

    std::vector<unsigned char> mode_image;
    bool first_level = true;

    while (true) {
        std::cout << "  Compressing " << current_w << "x" << current_h << "..." << std::endl;

        std::vector<char> compressed_level;

        if (first_level) {
            compress_image(current_image_data, current_w, current_h, compressed_level, &mode_image);
            first_level = false;
        } else {
            compress_image(current_image_data, current_w, current_h, compressed_level);
        }

        mip_data_levels.push_back(compressed_level);
        mip_sizes.push_back(compressed_level.size());
        total_compressed_size += compressed_level.size();

        if (current_w == 1 && current_h == 1) {
            break;
        }

        int next_w = std::max(1, current_w / 2);
        int next_h = std::max(1, current_h / 2);

        unsigned char* next_image_data = new unsigned char[next_w * next_h * 4];
        stbir_resize_uint8(current_image_data, current_w, current_h, 0,
            next_image_data, next_w, next_h, 0, 4);

        if (prev_image_data) {
            delete[] prev_image_data;
        }
        prev_image_data = current_image_data;
        current_image_data = next_image_data;

        current_w = next_w;
        current_h = next_h;
    }

    std::cout << "Done with mipmaps" << std::endl;

    if (prev_image_data) {
        delete[] prev_image_data;
    }
    delete[] current_image_data;
    //stbi_image_free(image_data);

    // Calculate dynamic metadata size: 9 base uint32_t fields + mip_offsets table
    // Metadata starts with 9 fields (indices 0-8), then mipmap offsets follow
    uint32_t metadata_u32_count = 9 + static_cast<uint32_t>(mip_sizes.size());
    uint32_t metadata_size = metadata_u32_count * sizeof(uint32_t);

    // Calculate mipmap offsets
    uint32_t current_offset = sizeof(TextureHeader) + metadata_size;
    for (size_t i = 0; i < mip_sizes.size(); ++i) {
        mip_offsets.push_back(current_offset);
        current_offset += mip_sizes[i];
    }


    // --- End of Mipmap Logic ---

    TextureHeader header;
    header.width = w;
    header.height = h;
    header.flags = 4;
    header.data_size = total_compressed_size;

    std::ofstream out_file(output_path, std::ios::binary);
    if (!out_file) {
        std::cerr << "Error: Cannot open output file: " << output_path << std::endl;
        return 1;
    }

    // Write the primary 24-byte header.
    out_file.write(reinterpret_cast<const char*>(&header), sizeof(header));

    // Create and populate the dynamic-sized metadata block.
    std::vector<char> metadata(metadata_size, 0);
    uint32_t* metadata_u32 = reinterpret_cast<uint32_t*>(metadata.data());

    if (transparent_image)
    {
        contains_alpha = 0x01;
        std::cout << "Setting the alpha bit to 0x01" << std::endl;
    }
    else
    {
        std::cout << "No alpha???" << std::endl;
    }

	metadata_u32[2] = contains_alpha;

    metadata_u32[3] = 0x20; // Compression flag
    metadata_u32[5] = 0x54455843; // "CXET"
    metadata_u32[7] = w;
    metadata_u32[8] = h;

    // Write the mipmap offset table, in reverse order (smallest mip's offset first)
    //TODO: check if this order is correct
    int table_index = 9;
    for (int i = mip_offsets.size() - 1; i >= 0; --i) {
        metadata_u32[table_index] = mip_offsets[i];
        table_index++;
    }

    //build the full file buffer of header, metadata to overwrite it using load_flags later
    std::vector<char> file_buffer;
    file_buffer.reserve(sizeof(TextureHeader) + metadata_size + total_compressed_size);

    const char* header_ptr = reinterpret_cast<const char*>(&header);
    file_buffer.insert(file_buffer.end(), header_ptr, header_ptr + sizeof(header));
    file_buffer.insert(file_buffer.end(), metadata.data(), metadata.data() + metadata.size());

    // Write the compressed data for each mip level, from largest to smallest.
    for (const auto& level_data : mip_data_levels) {
        file_buffer.insert(file_buffer.end(), level_data.data(), level_data.data() + level_data.size());
    }

    //std::string base_name;
    size_t dot_pos = input_path.find_last_of(".");
    if (dot_pos != std::string::npos) {
        base_name = input_path.substr(0, dot_pos);
    }
    else {
        base_name = input_path;
    }

    std::ifstream ini_check(ini_path);
    if (ini_check.good()) {
        std::cout << "Loading flags from " << ini_path << std::endl;
        load_flags(file_buffer.data(), ini_path);
    }

    if (!mode_image.empty()) {
        std::string modes_path = base_name + "_modes.png";
        if ((int)mode_image.size() == w * h * 4) {
            stbi_write_png(modes_path.c_str(), w, h, 4, mode_image.data(), w * 4);
            std::cout << "Saved debug mode map to " << modes_path << std::endl;
        }
    }

    out_file.seekp(0);
    out_file.write(file_buffer.data(), file_buffer.size());

    std::cout << "Successfully compressed and saved texture to " << output_path << std::endl;

    return 0;
}

int main(int argc, char* argv[]) {
	std::cout << "Mdk2VqTex (c) 2025 TurboSosiska a.k.a. Rusya" << std::endl;
	std::cout << "ubijca16@gmail.com" << std::endl << std::endl;
	
	bool do_flip = false;
    bool decompress_all_mips = false;

    if (argc < 3 || argc > 5) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file> [flip] [decompress_all_mips]" << std::endl;
		std::cerr << "input_file: can be either TEX file or PNG file with the texture that you want to process" << std::endl;
        std::cerr << "output_file: can be either TEX file or PNG file with the output texture" << std::endl;
        std::cerr << "flip: true/false - whether to perform a flip on Y axis" << std::endl;
        std::cerr << "decompress_all_mips: true/false - whether to decompress all the mipmaps" << std::endl;

        std::cerr << std::endl << "Examples:" << std::endl;
        std::cerr << "Mdk2VqTex texture.tex output.png true - this will convert texture.tex into output.png and flip it on Y axis" << std::endl;
        std::cerr << "Mdk2VqTex input.png output.tex true - this will convert input.png into texture.tex and flip it on Y axis" << std::endl;

        return 1;
    }
    std::string input_path = argv[1];
    std::string output_path = argv[2];

    std::string input_path_lower = input_path;
    std::transform(input_path_lower.begin(), input_path_lower.end(), input_path_lower.begin(),
        [](unsigned char c) { return std::tolower(c); });

    if (argc >= 4)
    {
        std::string flip_arg = argv[3];
        std::transform(flip_arg.begin(), flip_arg.end(), flip_arg.begin(),
            [](unsigned char c) { return std::tolower(c); });

        do_flip = flip_arg == "true";
    }

    if (argc == 5)
    {
        std::string mip_arg = argv[4];
        std::transform(mip_arg.begin(), mip_arg.end(), mip_arg.begin(),
            [](unsigned char c) { return std::tolower(c); });

        std::cout << "mip_arg: " << mip_arg << std::endl;

        decompress_all_mips = mip_arg == "true";
    }

    if (endsWith(input_path_lower, "png"))
    {
        return compress(input_path, output_path, do_flip, decompress_all_mips);
    }
    else
    {
        return decompress(input_path, output_path, do_flip, decompress_all_mips);
    }

    return 0;
}