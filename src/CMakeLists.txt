# executable target
set(TARGET qMandelbrot)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB SOURCE_FILES CONFIGURE_DEPENDS *.cpp *.h *.ui *.qrc LIST_DIRECTORIES false)
set(app_icon_macos "${CMAKE_CURRENT_SOURCE_DIR}/res/icons/Mandelbrot.icns")
set_source_files_properties(${app_icon_macos} PROPERTIES  MACOSX_PACKAGE_LOCATION "Resources")
set(app_icon_resource_windows "${CMAKE_CURRENT_SOURCE_DIR}/Mandelbrot.rc")

qt_add_executable(${TARGET} ${SOURCE_FILES} ${app_icon_macos} ${app_icon_resource_windows})

if(APPLE)
    # Build a macOS .app bundle with proper metadata and icon.
    set_target_properties(${TARGET} PROPERTIES
        MACOSX_BUNDLE TRUE
        MACOSX_BUNDLE_BUNDLE_NAME "qMandelbrot"
        MACOSX_BUNDLE_GUI_IDENTIFIER "com.ericgur.qMandelbrot"
        MACOSX_BUNDLE_BUNDLE_VERSION "1.0"
        MACOSX_BUNDLE_SHORT_VERSION_STRING "1.0"
        MACOSX_BUNDLE_ICON_FILE "Mandelbrot.icns"
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

install(TARGETS ${TARGET}
    BUNDLE  DESTINATION ${CMAKE_INSTALL_BINDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Deploy Qt runtime during cmake --install: on Windows copy DLLs/plugins next to
# the executable; on macOS macdeployqt embeds the Qt frameworks and plugins inside
# the .app bundle so it is self-contained.
if(WIN32 OR APPLE)
    message(STATUS "Configuring Qt deployment script...")
    qt_generate_deploy_app_script(
        TARGET ${TARGET}
        OUTPUT_SCRIPT deploy_script
        NO_TRANSLATIONS
        NO_UNSUPPORTED_PLATFORM_ERROR
    )
    if(WIN32)
        # CMAKE_INSTALL_PREFIX already points at <source>/bin, so Qt's deploy step should
        # treat the prefix itself as the bin dir instead of appending another "bin".
        install(CODE "set(QT_DEPLOY_BIN_DIR \".\")")
    endif()
    install(SCRIPT ${deploy_script})
endif()
