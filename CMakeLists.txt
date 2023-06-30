cmake_minimum_required(VERSION 2.8.12)

project(Ratatosk)

# To enable a larger default k-mer size, replace MAX_KMER_SIZE with a larger multiple of 32: actual maximum k-mer size will be MAX_KMER_SIZE-1.
SET(MAX_KMER_SIZE "64" CACHE STRING "MAX_KMER_SIZE")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE MATCHES Debug)
	message("Debug mode")
	add_compile_options(-g)
else(CMAKE_BUILD_TYPE MATCHES Debug)
	if(CMAKE_BUILD_TYPE MATCHES Profile)
		message("Profile mode")
		add_compile_options(-pg)
		set(CMAKE_SHARED_LINKER_FLAGS "-pg")
		set(CMAKE_EXE_LINKER_FLAGS "-pg")
	else(CMAKE_BUILD_TYPE MATCHES Profile)
		message("Release mode")
		add_compile_options(-O3)
	endif(CMAKE_BUILD_TYPE MATCHES Profile)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

add_subdirectory(src)
