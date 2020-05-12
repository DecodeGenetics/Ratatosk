file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp *.hh *.tcc)

# Ratatosk
add_definitions(-DMAX_KMER_SIZE=${MAX_KMER_SIZE})
add_executable(Ratatosk ${sources})

target_compile_options(Ratatosk PUBLIC "-DMAX_KMER_SIZE=${MAX_KMER_SIZE}")

# Add Bifrost as a dependency
target_link_libraries(Ratatosk bifrost)

# Add Pthread as a dependency
find_package(Threads REQUIRED)
target_link_libraries(Ratatosk ${CMAKE_THREAD_LIBS_INIT})

# Add zlib as a dependency
find_package(ZLIB REQUIRED)
target_link_libraries(Ratatosk ${ZLIB_LIBRARIES})

if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Zlib not found. Required for to output files")
endif(ZLIB_FOUND)

install(TARGETS Ratatosk DESTINATION bin)
