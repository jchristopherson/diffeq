# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(DIFFEQ_SOURCES
    ${dir}/diffeq.f90
    ${dir}/diffeq_base.f90
    ${dir}/diffeq_errors.f90
    ${dir}/diffeq_runge_kutta.f90
    ${dir}/diffeq_dprk45_constants.f90
    ${dir}/diffeq_bsrk32_constants.f90
    ${dir}/diffeq_rk853_constants.f90
    ${dir}/diffeq_rosenbrock_constants.f90
    ${dir}/diffeq_implicit_runge_kutta.f90
    ${dir}/VODE/dvode.f
    ${dir}/diffeq_vode.f90
)

# Dependencies
include(FetchContent)
find_package(ferror QUIET)
find_package(linalg QUIET)

if (NOT ferror_FOUND)
    message(STATUS "FERROR could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        ferror
        GIT_REPOSITORY "https://github.com/jchristopherson/ferror"
    )
    FetchContent_MakeAvailable(ferror)
    set(FERROR_LIBRARIES ferror)
    set(BUILD_FERROR TRUE)
else()
    set(FERROR_LIBRARIES ferror::ferror)
    set(BUILD_FERROR FALSE)
endif()

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
    ${FERROR_LIBRARIES}
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

if (${BUILD_FERROR} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${FERROR_LIBRARIES})
endif()

if (${BUILD_LINALG} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${LINALG_LIBRARIES})
endif()
