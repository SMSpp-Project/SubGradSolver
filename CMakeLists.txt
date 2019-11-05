# --------------------------------------------------------------------------- #
#    Main CMake file for SubGradSolver                                        #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.8...3.15)
# TODO: set proper version
project(SubGradSolver VERSION 1.0.0 LANGUAGES CXX)

# ----- Settings ------------------------------------------------------------ #
include(CTest)

# ----- Requirements -------------------------------------------------------- #
find_package(SMSpp REQUIRED)

# ----- Library ------------------------------------------------------------- #
add_library(SubGradSolver)
target_compile_features(SubGradSolver PUBLIC cxx_std_17)

target_sources(SubGradSolver PRIVATE
               Deflection.cpp
               LuckFi.cpp
               Stepsize.cpp
               SubGradSolver.cpp)

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake), that must be different
# for local build and installation.
target_include_directories(
        SubGradSolver PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/SMS++>)

target_link_libraries(SubGradSolver PUBLIC SMS++::SMSpp)

# Alias for testing
add_library(SMS++::SubGradSolver ALIAS SubGradSolver)

# ----- Install instructions ------------------------------------------------ #
# For installing the library and all its CMake modules on the system.
# Not required for local builds (see below).

include(GNUInstallDirs)

# Install library
install(TARGETS SubGradSolver
        EXPORT SubGradSolver-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install headers
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/Deflection.h
        ${CMAKE_CURRENT_SOURCE_DIR}/LuckFi.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Stepsize.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SubGradSolver.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Install target import file
install(EXPORT SubGradSolver-targets
        FILE SubGradSolverTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SubGradSolver)

# Generate package version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/SubGradSolverConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate package configuration file
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/SubGradSolverConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/SubGradSolverConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SubGradSolver)

# Install package version and configuration files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/SubGradSolverConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/SubGradSolverConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SubGradSolver)

# ----- Export from the build tree ------------------------------------------ #
# This is needed so other projects can find the library in a local build.
# It does not affect installation.

if (CMAKE_EXPORT_PACKAGE_REGISTRY)

    export(EXPORT SubGradSolver-targets
           FILE ${CMAKE_CURRENT_BINARY_DIR}/SubGradSolverTargets.cmake
           NAMESPACE SMS++::)

    # Register package in user's package registry (~/.cmake/registry)
    export(PACKAGE SubGradSolver)
endif ()

# ----- Subdirectories ------------------------------------------------------ #
if (BUILD_TESTING)
    add_subdirectory(test)
endif()
