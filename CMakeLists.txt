# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ SubGradSolver                                  #
#                                                                             #
#    This file allows one to build the library using CMake.                   #
#    To do so, you can use the following commands:                            #
#                                                                             #
#        $ cmake -S <source-path> -B <build-path>                             #
#        $ cmake --build <build-path>                                         #
#                                                                             #
#    The following command also installs the library in the system:           #
#                                                                             #
#        $ cmake --build <build-path> --target install                        #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.11)
project(SubGradSolver VERSION 0.2.0 LANGUAGES CXX)

# ----- Settings ------------------------------------------------------------ #
# This creates an ENABLE_TESTING option (default: ON) and enables the testing
# using the ctest executable.
# See: https://cmake.org/cmake/help/latest/manual/ctest.1.html
include(CTest)

# ----- Requirements -------------------------------------------------------- #
find_package(SMSpp REQUIRED)

# ----- Configuration header ------------------------------------------------ #
# This will generate a *Config.h header in the build directory.
configure_file(cmake/SubGradSolverConfig.h.in SubGradSolverConfig.h)

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(SubGradSolver)
target_compile_features(SubGradSolver PUBLIC cxx_std_17)

# When adding source files with target_sources(), PRIVATE means that the files
# should only be added to this library, whereas PUBLIC means they should be
# added to this library and to any target that links to it.
# INTERFACE can be used for sources that should not be added to this library
# but should be added to anything that links to it.
target_sources(SubGradSolver PRIVATE
               Deflection.cpp
               LuckFi.cpp
               Stepsize.cpp
               SubGradSolver.cpp)

# When using target_include_directories(), PUBLIC means that any targets
# that link to this target also need that include directory.
# Other options are PRIVATE (only affect the current target, not dependencies),
# and INTERFACE (only needed for dependencies).

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake).
# This means that if a target finds this library in its build directory
# will look into the BUILD_INTERFACE path for its headers, if it finds it
# installed in the system will look into the INSTALL_INTERFACE path.
target_include_directories(
        SubGradSolver PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/SMS++>)

# When linking other targets to the library with target_link_libraries(),
# PRIVATE means that the libraries will be linked only to this library,
# PUBLIC means they will be linked also to the targets that depend on this
# library, INTERFACE means they will be linked only to the targets that depend
# on this library.
target_link_libraries(SubGradSolver PUBLIC SMS++::SMSpp)

# This alias is defined so that executables in this same project can use
# the library with this notation.
add_library(SMS++::SubGradSolver ALIAS SubGradSolver)

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS SubGradSolver
        EXPORT SubGradSolver-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install the headers
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/Deflection.h
        ${CMAKE_CURRENT_SOURCE_DIR}/LuckFi.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Stepsize.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SubGradSolver.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/SubGradSolverConfig.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT SubGradSolver-targets
        FILE SubGradSolverTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SubGradSolver)

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/SubGradSolverConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate the package configuration file, that allows other
# CMake projects to find the library with find_package().
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/SubGradSolverConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/SubGradSolverConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SubGradSolver)

# Install the package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/SubGradSolverConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/SubGradSolverConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SubGradSolver)

# ----- Export from the build tree ------------------------------------------ #
# The following commands are needed so other projects can find the library
# in a local build directory. They use CMake's User Package Registry
# (located in ~/.cmake/registry) to do so.
# See: https://cmake.org/cmake/help/latest/manual/cmake-packages.7.html
# They don't affect the installation.

# CMAKE_EXPORT_PACKAGE_REGISTRY is OFF by default.
# It may be ON if the project is included by the SMS++ umbrella project.
if (CMAKE_EXPORT_PACKAGE_REGISTRY)

    # Generate the target import file, that allows other
    # CMake projects to import the target.
    export(EXPORT SubGradSolver-targets
           FILE ${CMAKE_CURRENT_BINARY_DIR}/SubGradSolverTargets.cmake
           NAMESPACE SMS++::)

    # Add the binary directory to RPATH, so executables that use
    # dynamic loading will look into it.
    set(SMSPP_BUILD_RPATH "${SMSPP_BUILD_RPATH};${CMAKE_CURRENT_BINARY_DIR}"
        PARENT_SCOPE)

    # Export the package in user's package registry
    export(PACKAGE SubGradSolver)
endif ()

# ----- Subdirectories ------------------------------------------------------ #
if (BUILD_TESTING)
    add_subdirectory(test)
endif ()
