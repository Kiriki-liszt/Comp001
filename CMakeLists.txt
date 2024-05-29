cmake_minimum_required(VERSION 3.14.0)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13 CACHE STRING "")

set(vst3sdk_SOURCE_DIR "C:/VST_SDK/vst3sdk")
if(NOT vst3sdk_SOURCE_DIR)
    message(FATAL_ERROR "Path to VST3 SDK is empty!")
endif()

project(Comp001
    # This is your plug-in version number. Change it here only.
    # Version number symbols usable in C++ can be found in
    # source/version.h and ${PROJECT_BINARY_DIR}/projectversion.h.
    VERSION 1.0.0.1
    DESCRIPTION "Comp001 VST 3 Plug-in"
)

set(SMTG_VSTGUI_ROOT "${vst3sdk_SOURCE_DIR}")

add_subdirectory(${vst3sdk_SOURCE_DIR} ${PROJECT_BINARY_DIR}/vst3sdk)
smtg_enable_vst3_sdk()

smtg_add_vst3plugin(Comp001
    source/version.h
    source/comp_cids.h
    source/comp_processor.h
    source/comp_processor.cpp
    source/comp_controller.h
    source/comp_controller.cpp
    source/comp_entry.cpp
    source/comp_dataexchange.h
)

#- VSTGUI Wanted ----
if(SMTG_ENABLE_VSTGUI_SUPPORT)
    target_sources(Comp001
        PRIVATE
            resource/comp_editor.uidesc
            resource/PretendardMono-Regular.ttf
    )
    target_link_libraries(Comp001
        PRIVATE
            vstgui_support
    )
    smtg_target_add_plugin_resources(Comp001
        RESOURCES
            "resource/comp_editor.uidesc"
            "resource/PretendardMono-Regular.ttf"
    )
endif(SMTG_ENABLE_VSTGUI_SUPPORT)
# -------------------

smtg_target_add_plugin_snapshots (Comp001
    RESOURCES
        resource/973A835973F55212BA5BA5476DAC162B_snapshot.png
        resource/973A835973F55212BA5BA5476DAC162B_snapshot_2.0x.png
)

target_link_libraries(Comp001
    PRIVATE
        sdk
)

smtg_target_configure_version_file(Comp001)

if(SMTG_MAC)
    smtg_target_set_bundle(Comp001
        BUNDLE_IDENTIFIER io.github.yg331.comp001
        COMPANY_NAME "yg331"
    )
    smtg_target_set_debug_executable(Comp001
        "/Applications/VST3PluginTestHost.app"
        "--pluginfolder;$(BUILT_PRODUCTS_DIR)"
    )
elseif(SMTG_WIN)
    target_sources(Comp001 PRIVATE 
        resource/win32resource.rc
    )
    if(MSVC)
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Comp001)

        smtg_target_set_debug_executable(Comp001
            "$(ProgramW6432)/Steinberg/VST3PluginTestHost/VST3PluginTestHost.exe"
            "--pluginfolder \"$(OutDir)/\""
        )
    endif()
endif(SMTG_MAC)
