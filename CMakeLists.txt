cmake_minimum_required(VERSION 3.5)

project (Mdk2VqTex)

set(SOURCES
    main.cpp
    compressor.cpp
    decompressor.cpp
)

set(HEADERS
    compressor.h
    decompressor.h
	
	stb_image.h
	stb_image_resize.h
	stb_image_write.h
)

add_executable(Mdk2VqTex ${SOURCES} ${HEADERS})

set_target_properties(Mdk2VqTex
    PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)