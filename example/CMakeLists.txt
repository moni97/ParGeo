add_executable(closestPair closestPair.cpp)
target_link_libraries(closestPair closestPairLib pargeoLib)

add_executable(graphGenerator graphGenerator.cpp)
target_link_libraries(graphGenerator
  pargeoLib
  delaunay2dLib
  emstLib
  spatialGraphLib)

add_executable(kNearestNeighbor kNearestNeighbor.cpp)
target_link_libraries(kNearestNeighbor pargeoLib)

add_executable(experimenting experimenting.cpp)
target_link_libraries(experimenting pargeoLib)

add_executable(readData readData.cpp)
target_link_libraries(readData pargeoLib)

add_executable(bichromaticClosestPair bichromaticClosestPair.cpp)
target_link_libraries(bichromaticClosestPair pargeoLib)

add_executable(rangeSearch rangeSearch.cpp)
target_link_libraries(rangeSearch pargeoLib)

add_executable(rangeSearchSampling rangeSearchSampling.cpp)
target_link_libraries(rangeSearchSampling pargeoLib)

add_executable(rangeEntropyTesting rangeEntropyTesting.cpp)
target_link_libraries(rangeEntropyTesting pargeoLib)

add_executable(rangeEntropyTesting2 rangeEntropyTesting2.cpp)
target_link_libraries(rangeEntropyTesting2 pargeoLib)

add_executable(getData getData.cpp)
target_link_libraries(getData pargeoLib)

add_executable(dynamicKdTreeKNN dynamicKdTreeKNN.cpp)
target_link_libraries(dynamicKdTreeKNN pargeoLib)

add_executable(batchKdTree batchKdTree.cpp)
target_link_libraries(batchKdTree pargeoLib)

add_executable(wellSeparatedPair wellSeparatedPair.cpp)
target_link_libraries(wellSeparatedPair pargeoLib)

add_executable(mortonSort mortonSort.cpp)
target_link_libraries(mortonSort pargeoLib)

add_executable(uniformDataGenerator uniformDataGenerator.cpp)
target_link_libraries(uniformDataGenerator pargeoLib)

add_executable(seedSpreader seedSpreaderGenerator.cpp)
target_link_libraries(seedSpreader
  datasetLib
  pargeoLib)

add_executable(convexHull convexHull.cpp)
target_link_libraries(convexHull pargeoLib
  divideHull2dLib
  pseudoHull3dLib
  )

add_executable(smallestEnclosingBall smallestEnclosingBall.cpp)
target_link_libraries(smallestEnclosingBall pargeoLib
  samplingSebLib
  # welzlSebLib
  # scanSebLib
  )
  