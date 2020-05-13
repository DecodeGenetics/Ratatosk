file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp *.hh *.tcc)

# Bifrost library
include(ExternalProject)

ExternalProject_Add(	
  bifrost_lib
  BUILD_IN_SOURCE 1
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/Bifrost
  PREFIX ${PROJECT_SOURCE_DIR}/Bifrost	
  CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DMAX_KMER_SIZE=${MAX_KMER_SIZE}
  UPDATE_COMMAND ""	
  INSTALL_COMMAND "")	

ExternalProject_Get_property(bifrost_lib INSTALL_DIR)

add_library(bifrost STATIC IMPORTED)	

set(BIFROST_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/Bifrost/src")
set(BIFROST_LIBRARY_DIR "${PROJECT_SOURCE_DIR}/Bifrost/src")

include_directories(SYSTEM ${BIFROST_INCLUDE_DIR})

set_property(TARGET bifrost PROPERTY IMPORTED_LOCATION "${BIFROST_LIBRARY_DIR}/libbifrost.a")

add_dependencies(bifrost bifrost_lib)

# Pthread library
find_package(Threads REQUIRED)

# Zlib library
find_package(ZLIB REQUIRED)

if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Zlib not found. Required for to output files")
endif(ZLIB_FOUND)

# Ratatosk
add_executable(Ratatosk ${sources})

# Ratatosk compilation
add_definitions(-DMAX_KMER_SIZE=${MAX_KMER_SIZE})
target_compile_options(Ratatosk PUBLIC "-DMAX_KMER_SIZE=${MAX_KMER_SIZE}")

# Ratatosk linking
target_link_libraries(Ratatosk ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(Ratatosk ${ZLIB_LIBRARIES})
target_link_libraries(Ratatosk bifrost)

install(TARGETS Ratatosk DESTINATION bin)
