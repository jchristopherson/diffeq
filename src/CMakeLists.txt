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
)
list(APPEND DIFFEQ_SOURCES ${DIFFEQ_C_SOURCES})


# Dependencies
include(FetchContent)
find_package(linalg 2.3.0 QUIET)

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
)
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

if (${BUILD_SHARED_LIBS})
    if (${BUILD_LINALG})
        install(TARGETS ${LINALG_LIBRARIES}
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
endif()
