set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${MED_INCLUDE_DIR}
  ${LIBMESH_INCLUDES} 
  ${MPI_INCLUDE_PATH}
  )

add_library(
  LIBMESHCXX SHARED 
  LibMeshUtil.cxx 
  LibMeshFunction.cxx
  LaplaceProblem.cxx 
  LinearElasticProblem.cxx 
  LIBMESH.cxx 
  _LibMeshProblem.cxx 
  ParallelMeshExtended.cxx
  EquationSystemsExtended.cxx
  Debug.cxx
  )

set(LIBMESHCXX_VERSION_MAJOR 1)
set(LIBMESHCXX_VERSION_MINOR 0)
set(LIBMESHCXX_VERSION_PATCH 1)
set(LIBMESHCXX_SOVERSION 
 "${LIBMESHCXX_VERSION_MAJOR}.${LIBMESHCXX_VERSION_MINOR}.${LIBMESHCXX_VERSION_PATCH}")

target_link_libraries(
  LIBMESHCXX
  ${MPI_LIBRARIES}
  ${MED_LIBRARIES}
  ${LIBMESH_LIBRARIES}
  )

set_target_properties(
  LIBMESHCXX
  PROPERTIES
  VERSION ${LIBMESHCXX_SOVERSION}
  SOVERSION ${LIBMESHCXX_VERSION_MAJOR})

install(
  TARGETS LIBMESHCXX
  DESTINATION lib)

install(
  FILES LIBMESH.hxx
  DESTINATION include)

IF (TEST_ENABLED)
  add_subdirectory(Test)
ENDIF (TEST_ENABLED)
