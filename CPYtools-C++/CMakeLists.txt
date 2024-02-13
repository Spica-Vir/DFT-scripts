# CMakeLists.txt
cmake_minimum_required(VERSION 3.11)
project(CRYSTALpytools_libs)

# C++
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXE_LINKER_FLAGS "-lstdc++ -lm")

# Python
find_package(Python3 COMPONENTS Development)

# Make libextfmt

set(SOURCES src/libextfmt.cpp src/libextfmt_python.cpp)
add_library(extfmt_lib SHARED ${SOURCES})

target_include_directories(extfmt_lib PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(extfmt_lib PRIVATE ${Python3_LIBRARIES})

set_target_properties(extfmt_lib PROPERTIES
	OUTPUT_NAME "extfmt"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib64"
)

install(TARGETS extfmt_lib DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/lib64")

