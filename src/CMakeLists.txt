# global properties
GET_PROPERTY(LOCAL_COMPILE_FLAGS GLOBAL PROPERTY MINORSEQ_COMPILE_FLAGS_GLOBAL)
GET_PROPERTY(LOCAL_LINK_FLAGS GLOBAL PROPERTY MINORSEQ_LINK_FLAGS_GLOBAL)

# get all header files for IDE support
file(GLOB_RECURSE UHU_HEADER "${UHU_IncludeDir}/*.h")
file(GLOB_RECURSE UHU_HIDDEN_HEADER "*.h")
# get sources for src/
file(GLOB UHU_CPP        "*.cpp")
# get sources for src/tools
file(GLOB UHU_TOOLS_CPP  "tools/*.cpp")

# ssw C library cannot be build with ${LOCAL_COMPILE_FLAGS} since it's C
add_library(ssw STATIC ${ssw_INCLUDE_DIRS}/ssw.c)

if (UHU_build_bin)
    function(create_exe exeName)
        add_executable(${exeName}
            ${UHU_SourceDir}/tools/main/${exeName}.cpp
            ${ssw_INCLUDE_DIRS}/ssw_cpp.cpp
        )

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
            ssw
        )

        target_include_directories(${exeName} PUBLIC
            ${UHU_IncludeDir}
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

    create_exe(zmw_stats)
    create_exe(lima)
endif()