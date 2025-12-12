set(diffeq_test_sources
    diffeq_test.f90
    diffeq_models.f90
    diffeq_jacobian_tests.f90
    diffeq_test_runge_kutta.f90
    diffeq_test_implicit_rk.f90
    diffeq_test_vode.f90
)

find_package(fortran_test_helper QUIET)
if (NOT fortran_test_helper_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        fortran_test_helper
        GIT_TAG "origin/main"
        GIT_REPOSITORY https://github.com/jchristopherson/fortran_test_helper
    )
    FetchContent_MakeAvailable(fortran_test_helper)
    set(FTH_LIBRARIES fortran_test_helper)
else()
    set(FTH_LIBRARIES fortran_test_helper::fortran_test_helper)
endif()

add_executable(diffeq_tests ${diffeq_test_sources})
target_link_libraries(
    diffeq_tests
    diffeq
    ${FTH_LIBRARIES}
)
add_test(
    NAME diffeq_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:diffeq_tests>
)

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET diffeq_tests
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_RUNTIME_DLLS:diffeq_tests>
        $<TARGET_FILE_DIR:diffeq_tests>
        COMMAND_EXPAND_LISTS
    )
endif()
