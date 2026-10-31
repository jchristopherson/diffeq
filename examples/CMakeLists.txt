# Get Dependencies
set(BUILD_TESTING OFF)
find_package(fplot QUIET)

if (fplot_FOUND)
    set(fplot_LIBRARY fplot::fplot)
else()
    include(FetchContent)
    FetchContent_Declare(
        fplot
        GIT_REPOSITORY "https://github.com/jchristopherson/fplot"
    )
    FetchContent_MakeAvailable(fplot)
    set(fplot_LIBRARY fplot)
endif()

# Shared model module used by the examples.
set(DIFFEQ_EXAMPLE_MODEL_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/model_mod_files)
add_library(
    diffeq_example_models
    STATIC
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
set_target_properties(
    diffeq_example_models PROPERTIES
    Fortran_MODULE_DIRECTORY ${DIFFEQ_EXAMPLE_MODEL_MOD_DIR}
)
target_include_directories(
    diffeq_example_models PUBLIC
    ${DIFFEQ_EXAMPLE_MODEL_MOD_DIR}
)

# RK45 Example
add_executable(
    rk45_example
    rk45_example.f90
)
target_link_libraries(rk45_example diffeq diffeq_example_models ${fplot_LIBRARY})

# RK853 Example
add_executable(
    rk853_example
    rk853_example.f90
)
target_link_libraries(rk853_example diffeq diffeq_example_models ${fplot_LIBRARY})

# Comparison Example
add_executable(
    comparison_example
    comparison_example.f90
)
target_link_libraries(comparison_example diffeq diffeq_example_models ${fplot_LIBRARY})

# Rosenbrock Example
add_executable(
    rosenbrock_example
    rosenbrock_example.f90
)
target_link_libraries(rosenbrock_example diffeq diffeq_example_models ${fplot_LIBRARY})

# Parametric Example
add_executable(
    parametric_example
    parametric_example.f90
)
target_link_libraries(parametric_example diffeq diffeq_example_models ${fplot_LIBRARY})

# DAE Example
add_executable(
    dae_example
    dae_example.f90
)
target_link_libraries(dae_example diffeq diffeq_example_models ${fplot_LIBRARY})

# Stiff Integrator Comparison Example
add_executable(
    stiff_integrator_example
    stiff_integrator_example.f90
)
target_link_libraries(stiff_integrator_example diffeq diffeq_example_models ${fplot_LIBRARY})

# Stabilized Pendulum Example
add_executable(
    stabilized_pendulum_example
    stabilized_pendulum_example.f90
)
target_link_libraries(stabilized_pendulum_example diffeq ${fplot_LIBRARY})

# ODE Results to File Example
add_executable(
    ode_results_to_file_example
    ode_results_to_file_example.f90
)
target_link_libraries(ode_results_to_file_example diffeq)

# Dense Output Example
add_executable(dense_output_example dense_output_example.f90)
target_link_libraries(dense_output_example diffeq diffeq_example_models ${fplot_LIBRARY})

# Event Detection Example
add_executable(event_detection_example event_detection_example.f90)
target_link_libraries(event_detection_example diffeq diffeq_example_models ${fplot_LIBRARY})

# AD Example
add_executable(ad_example ad_example.f90)
target_link_libraries(ad_example diffeq diffeq_example_models ${fplot_LIBRARY})

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET stiff_integrator_example
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:stiff_integrator_example> 
        $<TARGET_FILE_DIR:stiff_integrator_example>
        COMMAND_EXPAND_LISTS
    )
endif()

# ------------------------------------------------------------------------------
# C API Examples

add_executable(c_comparison_example c_comparison_example.c)
target_link_libraries(c_comparison_example diffeq)
target_include_directories(c_comparison_example PUBLIC ${CMAKE_SOURCE_DIR}/src/C)

add_executable(c_stabilized_pendulum_example c_stabilized_pendulum_example.c c_plot_helper.f90)
target_link_libraries(c_stabilized_pendulum_example diffeq ${fplot_LIBRARY})
target_include_directories(c_stabilized_pendulum_example PUBLIC ${CMAKE_SOURCE_DIR}/src/C ${CMAKE_SOURCE_DIR}/examples)

add_executable(c_vanderpol_example c_vanderpol_example.c c_plot_helper.f90)
target_link_libraries(c_vanderpol_example diffeq ${fplot_LIBRARY})
target_include_directories(c_vanderpol_example PUBLIC ${CMAKE_SOURCE_DIR}/src/C ${CMAKE_SOURCE_DIR}/examples)

add_executable(c_dense_output_example c_dense_output_example.c c_plot_helper.f90)
target_link_libraries(c_dense_output_example diffeq ${fplot_LIBRARY})
target_include_directories(c_dense_output_example PUBLIC ${CMAKE_SOURCE_DIR}/src/C ${CMAKE_SOURCE_DIR}/examples)

add_executable(c_event_detection_example c_event_detection_example.c c_plot_helper.f90)
target_link_libraries(c_event_detection_example diffeq ${fplot_LIBRARY})
target_include_directories(c_event_detection_example PUBLIC ${CMAKE_SOURCE_DIR}/src/C ${CMAKE_SOURCE_DIR}/examples)
