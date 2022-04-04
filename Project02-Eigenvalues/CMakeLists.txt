cmake_minimum_required(VERSION 3.12...3.22)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

set(NAME "test" CACHE STRING "The name of the project, `name.cpp` must be the main file.")

project(CompPhy_${NAME} LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(${CMAKE_BINARY_DIR}/conan_paths.cmake)
find_package(Eigen3)
find_package(yaml-cpp)

add_library(utils INTERFACE)
target_include_directories(utils INTERFACE include/utils)

add_executable(${NAME} ${NAME}.cpp)
target_link_libraries(${NAME} PUBLIC utils)
target_compile_features(${NAME} PUBLIC cxx_std_20)
target_link_libraries(${NAME} PUBLIC Eigen3::Eigen PUBLIC yaml-cpp)