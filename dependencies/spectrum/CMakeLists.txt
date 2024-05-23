include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    spectrum
    GIT_REPOSITORY "https://github.com/jchristopherson/spectrum"
    GIT_TAG main
)

FetchContent_MakeAvailable(spectrum)

set(spectrum_INCLUDE_DIR ${spectrum_BINARY_DIR}/include)
set(spectrum_INCLUDE_DIR ${spectrum_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(spectrum_LIBRARY spectrum)
set(spectrum_LIBRARY ${spectrum_LIBRARY} PARENT_SCOPE)

