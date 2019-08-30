cmake_minimum_required(VERSION 3.5.2)

option(test "Build all tests." OFF) # Makes boolean 'test' available.
SET(CMAKE_SKIP_BUILD_RPATH  TRUE)

set (CMAKE_CXX_STANDARD 11)
SET(CMAKE_C_COMPILER gcc ) 
SET(CMAKE_CXX_COMPILER g++ ) 
SET(CMAKE_CXX_FLAGS -pg)
SET(CMAKE_C_FLAGS -pg)

FIND_PATH(CONFIG++_INCLUDE_DIR libconfig.h++ /usr/include /usr/local/include)
FIND_LIBRARY(CONFIG++_LIBRARY NAMES config++ PATH /usr/lib /usr/local/lib) 
include_directories(${CONFIG++_INCLUDE_DIR})

FIND_PATH(MPI_INCLUDE_DIR mpi.h /usr/include /usr/local/include)
FIND_LIBRARY(MPI_LIBRARY NAMES mpi PATH /usr/lib /usr/local/lib) 
include_directories(${MPI_INCLUDE_DIR})

FIND_LIBRARY(PARDISO_LIBRARY NAMES pardiso PATH /usr/lib /usr/local/lib) 

FIND_PATH(BLAS_INCLUDE_DIR blas.h /usr/include /usr/local/include)
FIND_LIBRARY(BLAS_LIBRARY NAMES blas PATH /usr/lib /usr/local/lib) 
include_directories(${BLAS_INCLUDE_DIR})

project(LSweeps)

include(FindOpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

################################
# Normal Libraries & Executables
################################
find_package (Eigen3 REQUIRED NO_MODULE)

set( INCLUDE_DIR ./include/ )
include_directories( ${INCLUDE_DIR} )

set( SOURCE_DIR ./src/ )
set( SOURCEFILES 
  ${SOURCE_DIR}/FiniteDifferences.cpp 
  ${SOURCE_DIR}/AssemblerHelmholtz.cpp 
  ${SOURCE_DIR}/AssemblerHelmholtz2DFD.cpp
  ${SOURCE_DIR}/LocalLinearSolverEigenSparseLU.cpp
  ${SOURCE_DIR}/LocalLinearSolverPardiso.cpp
  ${SOURCE_DIR}/GlobalIterativeSolver.cpp
  ${SOURCE_DIR}/GlobalIterativeSolverGMRES.cpp
  ${SOURCE_DIR}/PML.cpp
  ${SOURCE_DIR}/PMLCubic2D.cpp
  ${SOURCE_DIR}/SubdomainStructured.cpp
  ${SOURCE_DIR}/DDCheckerboard.cpp
  ${SOURCE_DIR}/utils.cpp
  ${SOURCE_DIR}/Output.cpp
  ${SOURCE_DIR}/OutputVTK.cpp
  ${SOURCE_DIR}/Problem.cpp
  ${SOURCE_DIR}/Input.cpp
  ${SOURCE_DIR}/RHS.cpp
  ${SOURCE_DIR}/RHSSrcPts.cpp
  ${SOURCE_DIR}/Function.cpp
  ${SOURCE_DIR}/FunctionCSV.cpp
  ${SOURCE_DIR}/FunctionConstant.cpp
  )

add_library( lsweeps_lib ${SOURCEFILES})
target_link_libraries(lsweeps_lib Eigen3::Eigen ${MPI_LIBRARY} ${CONFIG++_LIBRARY} ${PARDISO_LIBRARY} ${BLAS_LIBRARY})

add_executable( lsweeps main.cpp )
target_link_libraries(lsweeps lsweeps_lib )

