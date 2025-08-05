cmake_minimum_required(VERSION 3.15)
project(DelSIEVE VERSION 1.4.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_XOPEN_SOURCE=500")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g -D_XOPEN_SOURCE=500")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")

# require the presence of some libraries
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

set(EXE_NAME DelSIEVE_simulator)

# set directories
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_library(rand SHARED ${SRC_DIR}/rand.cpp)
target_compile_options(rand PRIVATE -fPIC)

add_library(mtxexp SHARED ${SRC_DIR}/matrix_exponential.c ${SRC_DIR}/r8lib.c ${SRC_DIR}/c8lib.c)
target_compile_options(mtxexp PRIVATE -fPIC)

add_executable(${EXE_NAME}
        ${SRC_DIR}/cellcoal.c
        ${SRC_DIR}/eigen.c
        ${SRC_DIR}/hashMap.c
        ${SRC_DIR}/signatures.c
        )
target_link_libraries(${EXE_NAME} PRIVATE Threads::Threads OpenMP::OpenMP_C rand mtxexp)
if (UNIX AND NOT APPLE)
    target_link_libraries(${EXE_NAME} PRIVATE -lm -ldl)
endif ()

