cmake_minimum_required(VERSION 3.25)
project(qMandelbrot)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Pre-compiled header settings
set(CMAKE_PCH_INSTANTIATE_TEMPLATES ON)

# MacOS specific settings
set(CMAKE_OSX_DEPLOYMENT_TARGET 26.0)

find_package(OpenMP REQUIRED)

# set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -flto -march=native -Wall -Wextra -Wpedantic -Wimplicit-fallthrough -Wformat -Werror=format-security -Wno-gnu-anonymous-struct -Wno-nested-anon-types")
if (OpenMP_CXX_FOUND)
    link_libraries(OpenMP::OpenMP_CXX)
else()
    message(WARNING "OpenMP not found!")
endif()

# QT settings
find_package(Qt6 REQUIRED COMPONENTS Gui Widgets Core)
set(CMAKE_AUTORCC TRUE)
qt_standard_project_setup()

set(CMAKE_AUTOMOC_MOC_OPTIONS "-bpch.h")

#
# MacOS specific info configurations
#
if (APPLE)
    # Build a self-contained macOS .app bundle. Qt's deploy script (macdeployqt)
    # embeds the required Qt frameworks and plugins inside the bundle.
    set(CMAKE_MACOSX_BUNDLE TRUE)
    set(CMAKE_MACOSX_RPATH ON)
endif()

# Install everything under <source>/bin: the .app bundle on macOS, or the exe plus
# Qt DLLs/plugins on Windows. Keeping BINDIR "." stops Qt's deploy script from
# scattering plugins/ and translations/ into the source root.
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_INSTALL_BINDIR ".")
set(CMAKE_INSTALL_LIBDIR ".")

add_subdirectory(src)

# Windows: install compiler runtime for LLVM builds
if (QT_QMAKE_TARGET_MKSPEC STREQUAL "win32-clang-g++")
    # Install various runtime DLLs from the compiler's bin directory
    # This is needed for LLVM builds on Windows since qt deploy tool may or may not copy them.
    get_filename_component(TOOL_PATH ${CMAKE_CXX_COMPILER} DIRECTORY)
    install(FILES ${TOOL_PATH}/libomp.dll    DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(FILES ${TOOL_PATH}/libc++.dll    DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(FILES ${TOOL_PATH}/libunwind.dll DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
