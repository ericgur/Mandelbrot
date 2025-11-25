cmake_minimum_required(VERSION 3.5)
project(MandelbrotQt)

set(CMAKE_CXX_STANDARD 17)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/qt)

file(GLOB QT_SOURCES src/qt/*.cpp)

add_executable(MandelbrotQt ${QT_SOURCES})

target_link_libraries(MandelbrotQt Qt5::Widgets)
