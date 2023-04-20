# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(DIFFEQ_SOURCES
    ${dir}/diffeq.f90
)
set(DIFFEQ_SOURCES ${DIFFEQ_SOURCES} PARENT_SCOPE)