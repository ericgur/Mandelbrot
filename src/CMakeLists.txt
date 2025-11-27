# executable target
set(TARGET qMandelbrot)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB SOURCE_FILES CONFIGURE_DEPENDS *.cpp *.h *.ui *.qrc LIST_DIRECTORIES false)
set(app_icon_macos "${CMAKE_CURRENT_SOURCE_DIR}/res/icons/Mandelbrot.icns")
set_source_files_properties(${app_icon_macos} PROPERTIES  MACOSX_PACKAGE_LOCATION "Resources")
set(app_icon_resource_windows "${CMAKE_CURRENT_SOURCE_DIR}/Mandelbrot.rc")

qt_add_executable(${TARGET} ${SOURCE_FILES} ${app_icon_macos} ${app_icon_resource_windows})

if(APPLE)
    # Create a macOS .app bundle based on the global variable; build a plain executable and let the top-level CMake
    # collect and deploy runtimes into the single dist directory.
    set_target_properties(${TARGET} PROPERTIES
        MACOSX_BUNDLE ${CMAKE_MACOSX_BUNDLE}
    )
elseif(WIN32)
    set_target_properties(${TARGET} PROPERTIES
        WIN32_EXECUTABLE ON
    )
endif()

# precompiled header settings
target_precompile_headers(${TARGET} PUBLIC pch.h)

target_link_libraries(${TARGET} PRIVATE Qt6::Core Qt6::Widgets Qt6::Gui)
target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Qt6Gui_INCLUDE_DIRS} ${Qt6Widgets_INCLUDE_DIRS})

install(TARGETS ${TARGET} DESTINATION ${CMAKE_INSTALL_BINDIR})

# If on macOS set output directories and rpaths so everything is collected in ${DIST_DIR}
if(APPLE)
    # set output directories so all build artifacts go to the same folder
    set_target_properties(${TARGET} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${DIST_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${DIST_DIR}"
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        # use @loader_path so the binary will look next to itself for dylibs placed in the same directory
        INSTALL_RPATH "@loader_path"
        BUILD_RPATH "@loader_path"
    )
endif()
