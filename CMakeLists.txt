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
set(CMAKE_OSX_DEPLOYMENT_TARGET 15.0)

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
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

#
# MacOS specific info configurations
#
if (APPLE)
    # Do not create a .app bundle; collect all outputs into a single directory
    set(CMAKE_MACOSX_BUNDLE FALSE)

    # Distribution directory where all binaries and shared libs will be collected
    set(DIST_DIR "${CMAKE_BINARY_DIR}/dist")
    file(MAKE_DIRECTORY "${DIST_DIR}")

    # Keep FRAMEWORK_DIR variable for compatibility, point at DIST_DIR
    set(FRAMEWORK_DIR "${DIST_DIR}")
    set(CMAKE_INSTALL_BINDIR ${DIST_DIR})
    set(CMAKE_INSTALL_LIBDIR ${DIST_DIR})

    # Ensure rpath handling is enabled and prefer loader-relative paths at runtime
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_MACOSX_RPATH ON)
    set(CMAKE_INSTALL_RPATH "@loader_path")
    set(QT6_ROOT_DIR "$ENV{QT_DIR}")
    set(QT6_LIB_DIR "${QT6_ROOT_DIR}/macos/lib")    
endif()

add_subdirectory(src)

# Windows: install compiler runtime for LLVM builds
if (QT_QMAKE_TARGET_MKSPEC STREQUAL "win32-clang-g++")
    # install libc++.dll
    get_filename_component(TOOL_PATH ${CMAKE_CXX_COMPILER} DIRECTORY)
    install(FILES ${TOOL_PATH}/libc++.dll    DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
    install(FILES ${TOOL_PATH}/libunwind.dll DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
endif()

# Create symbolic links for Qt runtime libraries in the distribution directory.
# This places symlinks to the actual Qt6 library files (Core/Widgets/Gui)
# inside `DIST_DIR` so packaging/packagers can find them easily.
if (APPLE)
    message(STATUS "Creating symlinks for Qt6 frameworks in ${DIST_DIR}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${DIST_DIR}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${QT6_LIB_DIR}/QtCore.framework" "${DIST_DIR}/QtCore.framework")
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${QT6_LIB_DIR}/QtWidgets.framework" "${DIST_DIR}/QtWidgets.framework")
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${QT6_LIB_DIR}/QtGui.framework" "${DIST_DIR}/QtGui.framework")
endif()
