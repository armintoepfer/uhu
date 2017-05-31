
include(CheckCXXCompilerFlag)

# shared CXX flags for all source code & tests
set(UHU_FLAGS "-std=c++14 -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable")

# static linking
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(UHU_LINK_FLAGS "${UHU_LINK_FLAGS} -static-libgcc -static-libstdc++")
ENDIF()

# NOTE: quash clang warnings w/ Boost
check_cxx_compiler_flag("-Wno-unused-local-typedefs" HAS_NO_UNUSED_LOCAL_TYPEDEFS)
if(HAS_NO_UNUSED_LOCAL_TYPEDEFS)
    set(UHU_FLAGS "${UHU_FLAGS} -Wno-unused-local-typedefs")
endif()

# Cannot use this until pbbam complies
# if (CMAKE_COMPILER_IS_GNUCXX)
#     set(UHU_FLAGS "${UHU_FLAGS} -Werror=suggest-override")
# endif()

# Coverage settings
if (UHU_inc_coverage)
    set(UHU_FLAGS "${UHU_FLAGS} -fprofile-arcs -ftest-coverage")
endif()

# Extra testing that will lead to longer compilation times!
if (SANITIZE)
    # AddressSanitizer is a fast memory error detector
    set(UHU_FLAGS "${UHU_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")

    # Clang Thread Safety Analysis is a C++ language extension which warns about
    # potential race conditions in code.
    set(UHU_FLAGS "${UHU_FLAGS} -Wthread-safety")

    # ThreadSanitizer is a tool that detects data races
    set(UHU_FLAGS "${UHU_FLAGS} -fsanitize=thread")

    # MemorySanitizer is a detector of uninitialized reads.
    set(UHU_FLAGS "${UHU_FLAGS} -fsanitize=memory")

    # UndefinedBehaviorSanitizer is a fast undefined behavior detector.
    set(UHU_FLAGS "${UHU_FLAGS} -fsanitize=undefined")
endif()

if (JULIET_INHOUSE_PERFORMANCE)
    set(UHU_FLAGS "${UHU_FLAGS} -DJULIET_INHOUSE_PERFORMANCE")
endif()

# shared CXX flags for src & tests
SET_PROPERTY(GLOBAL PROPERTY MINORSEQ_COMPILE_FLAGS_GLOBAL ${UHU_FLAGS})
SET_PROPERTY(GLOBAL PROPERTY MINORSEQ_LINK_FLAGS_GLOBAL ${UHU_LINK_FLAGS})