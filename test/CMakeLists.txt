set(diffeq_test_sources
    diffeq_test.f90
    diffeq_models.f90
    diffeq_jacobian_tests.f90
    diffeq_test_runge_kutta.f90
    diffeq_test_implicit_rk.f90
    diffeq_test_bdf.f90
    diffeq_test_pece.f90
)

find_package(fortran_test_helper QUIET)
if (NOT fortran_test_helper_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        fortran_test_helper
        GIT_TAG "origin/main"
        GIT_REPOSITORY https://github.com/jchristopherson/fortran_test_helper
    )
    set(_linalg_build_testing ${BUILD_TESTING})
    set(BUILD_TESTING OFF)
    FetchContent_MakeAvailable(fortran_test_helper)
    set(BUILD_TESTING ${_linalg_build_testing})
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

add_executable(diffeq_error_tests diffeq_error_tests.f90)
target_link_libraries(diffeq_error_tests diffeq)

add_executable(diffeq_dense_output_tests diffeq_dense_output_tests.f90)
target_link_libraries(diffeq_dense_output_tests diffeq)
add_test(NAME diffeq_dense_output_tests COMMAND $<TARGET_FILE:diffeq_dense_output_tests>)
foreach(family singlestep multistep)
    foreach(error_case unavailable interval size nan infinite before after)
        add_test(
            NAME diffeq_dense_${family}_${error_case}
            COMMAND $<TARGET_FILE:diffeq_dense_output_tests> ${family}_${error_case}
        )
        set_tests_properties(diffeq_dense_${family}_${error_case} PROPERTIES WILL_FAIL TRUE)
    endforeach()

endforeach()

add_executable(diffeq_event_tests diffeq_event_tests.f90)
target_link_libraries(diffeq_event_tests diffeq)
add_test(NAME diffeq_event_tests COMMAND $<TARGET_FILE:diffeq_event_tests>)
foreach(error_case neqn equal_bounds y2_size y_size bracket_size bracket_index)
    add_test(
        NAME diffeq_event_error_${error_case}
        COMMAND $<TARGET_FILE:diffeq_event_tests> ${error_case}
    )
    set_tests_properties(diffeq_event_error_${error_case} PROPERTIES WILL_FAIL TRUE)
endforeach()

foreach(error_case missing_ode short_grid step_limit tolerance_size
    bdf_step_exhaustion adams_step_exhaustion
    kc_newton_step_limit kc_newton_tolerance
    ggl_dimensions ggl_gains ggl_state_size ggl_rhs_size ggl_matrix_size
    ggl_missing_callback ggl_tolerance ggl_singular_mass ggl_singular_constraints)
    add_test(
        NAME diffeq_error_${error_case}
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND $<TARGET_FILE:diffeq_error_tests> ${error_case}
    )
    set_tests_properties(diffeq_error_${error_case} PROPERTIES WILL_FAIL TRUE)
endforeach()

add_executable(diffeq_container_tests diffeq_container_tests.f90)
target_link_libraries(diffeq_container_tests diffeq)
add_test(
    NAME diffeq_container_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:diffeq_container_tests>
)

foreach(error_case
    init_zero_states init_negative_states init_zero_capacity init_negative_capacity
    capacity_uninitialized capacity_zero capacity_negative push_size
    point_uninitialized point_low point_high point_empty point_cleared
    selected_uninitialized selected_low selected_high selected_empty selected_cleared
    whole_uninitialized points_uninitialized
    vector_uninitialized vector_low vector_high vector_cleared
    matrix_uninitialized
)
    add_test(
        NAME diffeq_container_error_${error_case}
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND $<TARGET_FILE:diffeq_container_tests> ${error_case}
    )
    set_tests_properties(diffeq_container_error_${error_case} PROPERTIES WILL_FAIL TRUE)
endforeach()

add_executable(diffeq_container_io_tests diffeq_container_io_tests.f90)
target_link_libraries(diffeq_container_io_tests diffeq)
add_test(
    NAME diffeq_container_io_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:diffeq_container_io_tests>
)

add_executable(diffeq_c_api_tests diffeq_c_api_tests.c)
target_include_directories(diffeq_c_api_tests PRIVATE ${PROJECT_SOURCE_DIR}/src/C)
target_link_libraries(diffeq_c_api_tests diffeq)
add_test(
    NAME diffeq_c_api_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:diffeq_c_api_tests>
)

add_executable(diffeq_c_integrator_tests diffeq_c_integrator_tests.c)
target_include_directories(diffeq_c_integrator_tests PRIVATE ${PROJECT_SOURCE_DIR}/src/C)
target_link_libraries(diffeq_c_integrator_tests diffeq)
foreach(integrator rk23 rk45 rk853 tsitouras adams rosenbrock kc4 kc5 bdf)
    add_test(
        NAME diffeq_c_${integrator}
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND $<TARGET_FILE:diffeq_c_integrator_tests> ${integrator}
    )
endforeach()

add_executable(diffeq_c_rosenbrock_tests diffeq_c_rosenbrock_tests.c)
target_include_directories(diffeq_c_rosenbrock_tests PRIVATE ${PROJECT_SOURCE_DIR}/src/C)
target_link_libraries(diffeq_c_rosenbrock_tests diffeq)
add_test(
    NAME diffeq_c_rosenbrock_callbacks
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:diffeq_c_rosenbrock_tests>
)

add_executable(diffeq_c_dae_integrator_tests diffeq_c_dae_integrator_tests.c)
target_include_directories(diffeq_c_dae_integrator_tests PRIVATE ${PROJECT_SOURCE_DIR}/src/C)
target_link_libraries(diffeq_c_dae_integrator_tests diffeq)
foreach(integrator dae_rosenbrock dae_kc4 dae_kc5 dae_bdf)
    add_test(
        NAME diffeq_c_${integrator}
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND $<TARGET_FILE:diffeq_c_dae_integrator_tests> ${integrator}
    )
endforeach()

foreach(error_case quote carriage_return line_feed uninitialized open_error)
    add_test(
        NAME diffeq_container_io_error_${error_case}
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND $<TARGET_FILE:diffeq_container_io_tests> ${error_case}
    )
    set_tests_properties(diffeq_container_io_error_${error_case} PROPERTIES WILL_FAIL TRUE)
endforeach()

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET diffeq_tests
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_RUNTIME_DLLS:diffeq_tests>
        $<TARGET_FILE_DIR:diffeq_tests>
        COMMAND_EXPAND_LISTS
    )
    add_custom_command(
        TARGET diffeq_event_tests
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_RUNTIME_DLLS:diffeq_event_tests>
        $<TARGET_FILE_DIR:diffeq_event_tests>
        COMMAND_EXPAND_LISTS
    )
endif()
