@PACKAGE_INIT@

include(CMakeFindDependencyMacro)
find_dependency(linalg QUIET)
find_dependency(OpenMP REQUIRED COMPONENTS Fortran)
find_dependency(autodiff 1.0 CONFIG REQUIRED)

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
endif()