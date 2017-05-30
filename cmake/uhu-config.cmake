# Config generation

find_git_sha1(MINORSEQ_GIT_SHA1)

file (STRINGS "${UHU_RootDir}/CHANGELOG.md" MINORSEQ_CHANGELOG)

configure_file(
    ${UHU_IncludeDir}/pacbio/Version.h.in
    ${CMAKE_BINARY_DIR}/generated/pacbio/Version.h
)
