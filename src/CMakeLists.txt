# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(DIFFEQ_SOURCES
    ${dir}/diffeq.f90
    ${dir}/diffeq_base.f90
    ${dir}/diffeq_errors.f90
    ${dir}/containers/diffeq_containers.f90
    ${dir}/containers/diffeq_container_io.f90
    ${dir}/singlestep/diffeq_singlestep.f90
    ${dir}/singlestep/explicit/diffeq_dprk45_constants.f90
    ${dir}/singlestep/explicit/diffeq_tsit45_constants.f90
    ${dir}/singlestep/explicit/diffeq_bsrk32_constants.f90
    ${dir}/singlestep/explicit/diffeq_rk853_constants.f90
    ${dir}/singlestep/explicit/diffeq_tsitouras.f90
    ${dir}/singlestep/explicit/diffeq_rk45.f90
    ${dir}/singlestep/explicit/diffeq_rk23.f90
    ${dir}/singlestep/explicit/diffeq_rk853.f90
    ${dir}/singlestep/implicit/diffeq_rosenbrock_constants.f90
    ${dir}/singlestep/implicit/diffeq_rosenbrock.f90
    ${dir}/singlestep/implicit/diffeq_kennedy_carpenter.f90
    ${dir}/multistep/diffeq_multistep.f90
    ${dir}/multistep/diffeq_bdf.f90
    ${dir}/multistep/diffeq_pece.f90
    ${dir}/dae/diffeq_index_reduction.f90
)

set(DIFFEQ_C_SOURCES
    ${dir}/C/diffeq_c_helpers.f90
    ${dir}/C/diffeq_c_errors.f90
    ${dir}/C/diffeq_c_containers.f90
    ${dir}/C/diffeq_c_integrators.f90
    ${dir}/C/diffeq_c_index_reduction.f90
    ${dir}/C/diffeq_c_function_pointers.f90
)
list(APPEND DIFFEQ_SOURCES ${DIFFEQ_C_SOURCES})


# Dependencies
include(FetchContent)
find_package(OpenMP REQUIRED COMPONENTS Fortran)
find_package(autodiff 1.0 CONFIG QUIET)
if (NOT autodiff_FOUND)
    message(STATUS "AUTODIFF could not be found. A reference version will be employed.")
    FetchContent_Declare(
        autodiff
        GIT_TAG "origin/main"
        GIT_REPOSITORY "https://github.com/jchristopherson/autodiff"
    )
    set(_diffeq_build_testing ${BUILD_TESTING})
    set(BUILD_TESTING OFF)
    set(AUTODIFF_ENABLE_OPENMP ON)
    FetchContent_MakeAvailable(autodiff)
    set(BUILD_TESTING ${_diffeq_build_testing})
    set(AUTODIFF_LIBRARIES autodiff)
    set(BUILD_AUTODIFF TRUE)
else()
    set(AUTODIFF_LIBRARIES autodiff::autodiff)
    set(BUILD_AUTODIFF FALSE)
endif()
find_package(linalg 3.0 QUIET)

if (NOT linalg_FOUND)
    message(STATUS "LINALG could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        linalg
        GIT_REPOSITORY "https://github.com/jchristopherson/linalg"
    )
    FetchContent_MakeAvailable(linalg)
    set(LINALG_LIBRARIES linalg)
    set(BUILD_LINALG TRUE)
else()
    set(LINALG_LIBRARIES linalg::linalg)
    set(BUILD_LINALG FALSE)
endif()

# Build
add_library(${PROJECT_NAME} ${DIFFEQ_SOURCES})
target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    ${LINALG_LIBRARIES}
    ${AUTODIFF_LIBRARIES}
)
if (OpenMP_Fortran_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran)
    target_compile_definitions(${PROJECT_NAME} PRIVATE USEOPENMP=1)

    if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-fopenmp>)
        target_link_options(${PROJECT_NAME} PUBLIC $<$<LINK_LANGUAGE:Fortran>:-fopenmp>)
    elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
        target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-qopenmp>)
        target_link_options(${PROJECT_NAME} PUBLIC $<$<LINK_LANGUAGE:Fortran>:-qopenmp>)
    elseif (MSVC)
        target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:/openmp>)
        target_link_options(${PROJECT_NAME} PUBLIC $<$<LINK_LANGUAGE:Fortran>:/openmp>)
    endif()
endif()
set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)

set(DIFFEQ_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files/)
if (NOT EXISTS "${DIFFEQ_MOD_DIR}")
    file(MAKE_DIRECTORY "${DIFFEQ_MOD_DIR}")
endif()

set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY ${DIFFEQ_MOD_DIR}
)
target_include_directories(
    ${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${DIFFEQ_MOD_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>
)

# Install
install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${DIFFEQ_MOD_DIR} DESTINATION ${CMAKE_INSTALL_MODULEDIR})
install(
    FILES ${dir}/C/diffeq.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

if (${BUILD_SHARED_LIBS})
    if (${BUILD_LINALG})
        install(TARGETS ${LINALG_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_AUTODIFF})
        install(TARGETS ${AUTODIFF_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
else()
    if (${BUILD_LINALG})
        install(TARGETS ${LINALG_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_AUTODIFF})
        install(TARGETS ${AUTODIFF_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
endif()
