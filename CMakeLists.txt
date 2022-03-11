project(mpm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

cmake_minimum_required(VERSION 3.12)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# General compile settings
IF (NOT CMAKE_BUILD_TYPE)
 #SET(CMAKE_BUILD_TYPE "Debug")
 SET(CMAKE_BUILD_TYPE "Release")
ENDIF (NOT CMAKE_BUILD_TYPE)

# GNU specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()

# Intel specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

# Clang specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
endif()

# Check if Ninja
set(USED_CMAKE_GENERATOR "${CMAKE_GENERATOR}" CACHE STRING "Expose CMAKE_GENERATOR" FORCE)
if (USED_CMAKE_GENERATOR MATCHES "Ninja")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
endif()

# CMake seems to have no way to enable/disable testing per subproject,
# so we provide an option similar to BUILD_TESTING, but just for MPM.
option(MPM_BUILD_TESTING "enable testing for mpm" ON)

# Halo exchange
option(HALO_EXCHANGE "Enable halo exchange" OFF)

# PETSC
option(USE_PETSC "Use PETSC solver library" OFF)

# CMake Modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Boost Archive
find_package(Boost REQUIRED COMPONENTS filesystem system)
include_directories(${BOOST_INCLUDE_DIRS})
link_libraries(${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# MPI
find_package(MPI)
if (MPI_FOUND)
  if(HALO_EXCHANGE)
    add_definitions(-DUSE_HALO_EXCHANGE)
  endif()
  add_definitions("-DUSE_MPI")
  include_directories(${MPI_CXX_INCLUDE_DIRS})
  link_libraries(${MPI_CXX_LIBRARIES})
endif()

# HDF5
ENABLE_LANGUAGE(C)
find_package(HDF5 COMPONENTS CXX HL)
if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
  link_libraries(${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES})
  add_definitions(${HDF5_DEFINITIONS})
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof")
endif()

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
  if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

# pthreads
find_package (Threads)
link_libraries (${CMAKE_THREAD_LIBS_INIT})

# MKL
find_package(MKL)
if (MKL_FOUND)
  if (NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    link_libraries(iomp5)
    include_directories(${MKL_INCLUDE_DIR})
    link_libraries(${MKL_LIBRARIES})
    add_definitions("-DUSE_MKL")
  endif()
endif()

# PETSc
if (USE_PETSC)
  find_package(PETSc COMPONENTS C)
endif()
if (PETSC_FOUND)
  include_directories(${PETSC_INC})
  link_libraries(${PETSC_LIB})
  add_definitions("-DUSE_PETSC")
endif()

# KaHIP
if (MPI_FOUND)
  if (NO_KAHIP)
    find_package(KaHIP)
  else()
    find_package(KaHIP REQUIRED)
  endif()
if (KAHIP_FOUND)
  add_definitions("-DUSE_GRAPH_PARTITIONING")
  include_directories(${KAHIP_INCLUDE_DIRS})
  link_libraries(${KAHIP_LIBRARIES})
endif()
endif()

# VTK
find_package(VTK)
if (VTK_FOUND)
  add_definitions("-DUSE_VTK")
  include_directories(${VTK_INCLUDE_DIRS})
  link_libraries(${VTK_LIBRARIES})
endif()

# Partio
find_package(Partio)
if (PARTIO_FOUND)
  add_definitions("-DUSE_PARTIO")
  include_directories(${PARTIO_INCLUDE_DIRS})
  link_libraries(${PARTIO_LIBRARIES})
endif()

# Include directories
include_directories(BEFORE
  ${mpm_SOURCE_DIR}/include/
  ${mpm_SOURCE_DIR}/include/cells/
  ${mpm_SOURCE_DIR}/include/contacts/
  ${mpm_SOURCE_DIR}/include/containers/
  ${mpm_SOURCE_DIR}/include/data_structures/
  ${mpm_SOURCE_DIR}/include/elements/
  ${mpm_SOURCE_DIR}/include/elements/2d
  ${mpm_SOURCE_DIR}/include/elements/3d
  ${mpm_SOURCE_DIR}/include/functions/
  ${mpm_SOURCE_DIR}/include/generators/
  ${mpm_SOURCE_DIR}/include/io/
  ${mpm_SOURCE_DIR}/include/linear_solvers/
  ${mpm_SOURCE_DIR}/include/linear_solvers/assemblers
  ${mpm_SOURCE_DIR}/include/linear_solvers/convergence_criteria
  ${mpm_SOURCE_DIR}/include/linear_solvers/linear_solvers
  ${mpm_SOURCE_DIR}/include/loads_bcs/
  ${mpm_SOURCE_DIR}/include/materials/
  ${mpm_SOURCE_DIR}/include/mesh/
  ${mpm_SOURCE_DIR}/include/nodes/
  ${mpm_SOURCE_DIR}/include/particles/
  ${mpm_SOURCE_DIR}/include/particles/pod_particles
  ${mpm_SOURCE_DIR}/include/solvers/
  ${mpm_SOURCE_DIR}/include/solvers/mpm_scheme/
  ${mpm_SOURCE_DIR}/include/utilities/
  ${mpm_SOURCE_DIR}/external/
  ${mpm_SOURCE_DIR}/tests/include/
)

# mpm executable
SET(mpm_src
  ${mpm_SOURCE_DIR}/src/affine_transform.cc
  ${mpm_SOURCE_DIR}/src/cell.cc
  ${mpm_SOURCE_DIR}/src/element.cc
  ${mpm_SOURCE_DIR}/src/functions/functions.cc
  ${mpm_SOURCE_DIR}/src/functions/linear_function.cc
  ${mpm_SOURCE_DIR}/src/functions/sin_function.cc
  ${mpm_SOURCE_DIR}/src/geometry.cc
  ${mpm_SOURCE_DIR}/src/pod_particle.cc
  ${mpm_SOURCE_DIR}/src/pod_particle_twophase.cc
  ${mpm_SOURCE_DIR}/src/io/io.cc
  ${mpm_SOURCE_DIR}/src/io/io_mesh.cc
  ${mpm_SOURCE_DIR}/src/io/logger.cc
  ${mpm_SOURCE_DIR}/src/io/partio_writer.cc
  ${mpm_SOURCE_DIR}/src/io/vtk_writer.cc
  ${mpm_SOURCE_DIR}/src/linear_solver.cc
  ${mpm_SOURCE_DIR}/src/material.cc
  ${mpm_SOURCE_DIR}/src/mpm.cc
  ${mpm_SOURCE_DIR}/src/nodal_properties.cc
  ${mpm_SOURCE_DIR}/src/node.cc
  ${mpm_SOURCE_DIR}/src/particle.cc
  ${mpm_SOURCE_DIR}/src/quadrature.cc
)
add_executable(mpm ${mpm_SOURCE_DIR}/src/main.cc ${mpm_src} ${mpm_vtk})

# Git revision
set(PRE_CONFIGURE_FILE "${CMAKE_SOURCE_DIR}/include/git.cc.in")
set(POST_CONFIGURE_FILE "${CMAKE_SOURCE_DIR}/include/git.cc")
include("${CMAKE_SOURCE_DIR}/cmake/git_watcher.cmake")
add_library(git STATIC ${POST_CONFIGURE_FILE})
target_include_directories(git PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_dependencies(git check_git)
target_link_libraries(mpm git)

# Unit test
if(MPM_BUILD_TESTING)
  SET(test_src
    ${mpm_SOURCE_DIR}/tests/test_main.cc
    ${mpm_SOURCE_DIR}/tests/materials/bound_surf_plasticity_tx.cc
  )
  add_executable(mpmtest ${mpm_src} ${test_src})
  add_test(NAME mpmtest COMMAND $<TARGET_FILE:mpmtest>)
  enable_testing()

endif()

# Coverage
find_package(codecov)
if(ENABLE_COVERAGE)
  add_executable(mpmtest_coverage ${mpm_src} ${test_src})
  add_coverage(mpmtest_coverage)
endif()
