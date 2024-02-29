include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    nonlin
    GIT_REPOSITORY "https://github.com/jchristopherson/nonlin"
)

FetchContent_MakeAvailable(nonlin)

set(nonlin_INCLUDE_DIR ${nonlin_BINARY_DIR}/include)
set(nonlin_INCLUDE_DIR ${nonlin_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(nonlin_LIBRARY nonlin)
set(nonlin_LIBRARY ${nonlin_LIBRARY} PARENT_SCOPE)
