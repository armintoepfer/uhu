# global properties
GET_PROPERTY(LOCAL_COMPILE_FLAGS GLOBAL PROPERTY MINORSEQ_COMPILE_FLAGS_GLOBAL)
GET_PROPERTY(LOCAL_LINK_FLAGS GLOBAL PROPERTY MINORSEQ_LINK_FLAGS_GLOBAL)

# get all header files for IDE support
file(GLOB_RECURSE MS_HEADER "${MS_IncludeDir}/*.h")
file(GLOB_RECURSE MS_HIDDEN_HEADER "*.h")
# get sources for src/
file(GLOB MS_CPP        "*.cpp")
# get sources for src/tools
file(GLOB MS_TOOLS_CPP  "tools/*.cpp")

if (MS_build_bin)
    function(create_exe exeName)
        add_executable(${exeName} ${MS_SourceDir}/tools/main/${exeName}.cpp)

        set_target_properties(${exeName} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
        )

        target_link_libraries(${exeName}
            ${CMAKE_THREAD_LIBS_INIT}
            ${CMAKE_DL_LIBS}
            ${ZLIB_LIBRARIES}
            ${HTSLIB_LIBRARIES}
            ${PacBioBAM_LIBRARIES}
            ${pbcopper_LIBRARIES}
        )

        target_include_directories(${exeName} PUBLIC
            ${MS_IncludeDir}
            ${CMAKE_BINARY_DIR}/generated
            ${Boost_INCLUDE_DIRS}
            ${ZLIB_INCLUDE_DIRS}
            ${HTSLIB_INCLUDE_DIRS}
            ${PacBioBAM_INCLUDE_DIRS}
            ${pbcopper_INCLUDE_DIRS}
            ${ssw_INCLUDE_DIRS}
        )

        install(TARGETS ${exeName} RUNTIME DESTINATION bin)
        set_target_properties(${exeName} PROPERTIES COMPILE_FLAGS ${LOCAL_COMPILE_FLAGS})
        if (LOCAL_LINK_FLAGS)
            set_target_properties(${exeName} PROPERTIES LINK_FLAGS ${LOCAL_LINK_FLAGS})
        endif()
    endfunction()

    create_exe(extract_cx)
endif()