include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    fftpack
    GIT_REPOSITORY "https://github.com/fortran-lang/fftpack"
    GIT_TAG main
)

FetchContent_MakeAvailable(fftpack)

set(fftpack_INCLUDE_DIR ${fftpack_BINARY_DIR}/include)
set(fftpack_INCLUDE_DIR ${fftpack_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(fftpack_LIBRARY fftpack)
set(fftpack_LIBRARY ${fftpack_LIBRARY} PARENT_SCOPE)

