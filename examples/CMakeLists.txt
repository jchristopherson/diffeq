# Get Dependencies
find_package(fplot QUIET)

if (fplot_FOUND)
    # The installed fplot package can export transitive imported targets
    # like ferror::ferror, collections::collections, geompack::geompack,
    # and fstring::fstring. Ensure those dependency targets are imported
    # before linking so they are not passed through as raw -l flags.
    find_package(ferror QUIET)
    find_package(collections QUIET)
    find_package(geompack QUIET)
    find_package(fstring QUIET)

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

# RK45 Example
add_executable(
    rk45_example
    rk45_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(rk45_example diffeq ${fplot_LIBRARY})

# RK853 Example
add_executable(
    rk853_example
    rk853_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(rk853_example diffeq ${fplot_LIBRARY})

# Comparison Example
add_executable(
    comparison_example
    comparison_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(comparison_example diffeq ${fplot_LIBRARY})

# Rosenbrock Example
add_executable(
    rosenbrock_example
    rosenbrock_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(rosenbrock_example diffeq ${fplot_LIBRARY})

# Parametric Example
add_executable(
    parametric_example
    parametric_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(parametric_example diffeq ${fplot_LIBRARY})

# Stiff Integrator Comparison Example
add_executable(
    stiff_integrator_example
    stiff_integrator_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(stiff_integrator_example diffeq ${fplot_LIBRARY})

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