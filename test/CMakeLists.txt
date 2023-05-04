include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

add_subdirectory(fortran_test_helper)

set(diffeq_test_sources
    diffeq_test.f90
    diffeq_models.f90
    diffeq_jacobian_tests.f90
    diffeq_rkfixed_tests.f90
    diffeq_expfixed_tests.f90
    diffeq_adamsfixed_tests.f90
    diffeq_dprk45_tests.f90
)

add_executable(diffeq_tests ${diffeq_test_sources})
link_library(diffeq_tests ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(diffeq_tests ${fortran_test_helper_LIBRARY} ${fortran_test_helper_INCLUDE_DIR})
add_test(
    NAME diffeq_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:diffeq_tests>
)