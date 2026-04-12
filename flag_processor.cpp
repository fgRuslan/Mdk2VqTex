#include "flag_processor.h"
#include <iostream>
#include "SimpleIni.h"

void dump_flags(char* texture, std::string output_file)
{
	CSimpleIniA ini;
	ini.SetUnicode();

	uint8_t mipmap = *((uint32_t*)texture + 5);
	uint8_t fixed = *((uint32_t*)texture + 6);
	uint8_t clamp = *((uint32_t*)texture + 7);
	uint8_t blendmode = *((uint32_t*)texture + 8);
	uint8_t compression = (*((uint32_t*)texture + 9) & 0x20) >> 5;

	ini.SetLongValue("flags", "mipmap", mipmap);
	ini.SetLongValue("flags", "fixed", fixed);
	ini.SetLongValue("flags", "clamp", clamp);
	ini.SetLongValue("flags", "blendmode", blendmode);
	ini.SetLongValue("flags", "VQ compression", compression);

	ini.SaveFile(output_file.c_str());
}

void load_flags(char* texture, std::string input_file)
{
	CSimpleIniA ini;
	ini.SetUnicode();
	ini.LoadFile(input_file.c_str());

	uint8_t mipmap = ini.GetLongValue("flags", "mipmap", 0);
	uint8_t fixed = ini.GetLongValue("flags", "fixed", 0);
	uint8_t clamp = ini.GetLongValue("flags", "clamp", 0);
	uint8_t blendmode = ini.GetLongValue("flags", "blendmode", 0);
	uint8_t compression = ini.GetLongValue("flags", "VQ compression", 0x20);

	*((uint32_t*)texture + 5) = mipmap;
	*((uint32_t*)texture + 6) = fixed;
	*((uint32_t*)texture + 7) = clamp;
	*((uint32_t*)texture + 8) = blendmode;

	uint32_t compressionWord = *((uint32_t*)texture + 9);
	compressionWord = (compressionWord & ~0x20) | (compression << 5);
	*((uint32_t*)texture + 9) = compressionWord;
}