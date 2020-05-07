file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp *.hh *.tcc)

# Ratatosk

add_definitions(-DMAX_KMER_SIZE=${MAX_KMER_SIZE})
add_executable(Ratatosk ${sources})

target_compile_options(Ratatosk PUBLIC "-DMAX_KMER_SIZE=${MAX_KMER_SIZE}")

# Add Bifrost as a dependency

include(ExternalProject)

ExternalProject_Add(
  bifrost
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/Bifrost
  CMAKE_ARGS "${CMAKE_ARGS};-DMAX_KMER_SIZE=${MAX_KMER_SIZE} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")

ExternalProject_Get_property(bifrost INSTALL_DIR)

set(bifrost_INCLUDE "${CMAKE_SOURCE_DIR}/Bifrost/src")
set(bifrost_LIB "${INSTALL_DIR}/src/bifrost-build/src")

add_dependencies(Ratatosk bifrost)

target_include_directories(Ratatosk PUBLIC "${bifrost_INCLUDE}")
target_link_libraries(Ratatosk "${bifrost_LIB}/libbifrost.so")

# Add Pthread as a dependcy

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
