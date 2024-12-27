cmake_minimum_required (VERSION 3.9)

project (GeoVec-VE
  VERSION 1.0.0
  LANGUAGES CXX
)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

find_package(PNG REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS system)
find_library(MAPNIK mapnik)
find_library(HIGIS NAMES libhigis.a libhigis.so PATHS ../../Algorithms/library/dep-lib NO_DEFAULT_PATH)

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS"${CMAKE_EXE_LINKER_FLAGS}${OpenMP_EXE_LINKER_FLAGS}")
endif()


add_definitions(-std=c++17 -Wall -O3 -g)

set(INCLUDE_DIR /usr/local/include ../../Algorithms/include ../../Algorithms/include/dep-include)
include_directories(${INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})


add_executable(exp1_fig8  
              ../../Algorithms/library/buildIndex.cpp
              ../../Algorithms/library/indexNode.cpp 
              ../../Algorithms/library/topoJudge.cpp 
              ../../Algorithms/library/pointRender.cpp
              ../../Algorithms/library/linestringRender.cpp
              ../../Algorithms/library/polygonRender.cpp
              ${PROJECT_SOURCE_DIR}/source-code/exp1_fig8.cpp
              )
target_link_libraries(exp1_fig8 ${HIGIS} ${PNG_LIBRARIES} ${Boost_LIBRARIES} ${MAPNIK} ${MPI_LIBRARIES} -lmpich -lpthread)


add_executable(exp1_fig9  
              ../../Algorithms/library/buildIndex.cpp
              ../../Algorithms/library/indexNode.cpp 
              ../../Algorithms/library/topoJudge.cpp 
              ../../Algorithms/library/pointRender.cpp
              ../../Algorithms/library/linestringRender.cpp
              ../../Algorithms/library/polygonRender.cpp
              ${PROJECT_SOURCE_DIR}/source-code/exp1_fig9.cpp
              )
target_link_libraries(exp1_fig9 ${HIGIS} ${PNG_LIBRARIES} ${Boost_LIBRARIES} ${MAPNIK} ${MPI_LIBRARIES} -lmpich -lpthread)


add_executable(exp1_fig10  
              ../../Algorithms/library/buildIndex.cpp
              ../../Algorithms/library/indexNode.cpp 
              ../../Algorithms/library/topoJudge.cpp 
              ../../Algorithms/library/linestringRender.cpp
              ../../Algorithms/library/polygonRender.cpp
              ${PROJECT_SOURCE_DIR}/source-code/exp1_fig10.cpp
              )
target_link_libraries(exp1_fig10 ${HIGIS} ${PNG_LIBRARIES} ${Boost_LIBRARIES} ${MAPNIK} ${MPI_LIBRARIES} -lmpich -lpthread)


add_executable(exp1_fig11  
              ../../Algorithms/library/buildIndex.cpp
              ../../Algorithms/library/indexNode.cpp 
              ../../Algorithms/library/topoJudge.cpp 
              ../../Algorithms/library/polygonRender.cpp
              ${PROJECT_SOURCE_DIR}/source-code/exp1_fig11.cpp
              )
target_link_libraries(exp1_fig11 ${HIGIS} ${PNG_LIBRARIES} ${Boost_LIBRARIES} ${MAPNIK} ${MPI_LIBRARIES} -lmpich -lpthread)


add_executable(exp1_fig12  
              ../../Algorithms/library/buildIndex.cpp
              ../../Algorithms/library/indexNode.cpp 
              ../../Algorithms/library/topoJudge.cpp 
              ../../Algorithms/library/pointRender.cpp
              ../../Algorithms/library/linestringRender.cpp
              ../../Algorithms/library/polygonRender.cpp
              ${PROJECT_SOURCE_DIR}/source-code/exp1_fig12.cpp
              )
target_link_libraries(exp1_fig12 ${HIGIS} ${PNG_LIBRARIES} ${Boost_LIBRARIES} ${MAPNIK} ${MPI_LIBRARIES} -lmpich -lpthread)


add_executable(exp2_fig13-14  
              ../../Algorithms/library/buildIndex.cpp
              ../../Algorithms/library/indexNode.cpp 
              ../../Algorithms/library/topoJudge.cpp 
              ../../Algorithms/library/pointRender.cpp
              ../../Algorithms/library/linestringRender.cpp
              ../../Algorithms/library/polygonRender.cpp
              ${PROJECT_SOURCE_DIR}/source-code/exp2_fig13-14.cpp
              )
target_link_libraries(exp2_fig13-14 ${HIGIS} ${PNG_LIBRARIES} ${Boost_LIBRARIES} ${MAPNIK} ${MPI_LIBRARIES} -lmpich -lpthread)


add_executable(exp3_fig15  
              ../../Algorithms/library/buildIndex.cpp
              ../../Algorithms/library/indexNode.cpp 
              ../../Algorithms/library/topoJudge.cpp 
              ../../Algorithms/library/pointRender.cpp
              ../../Algorithms/library/linestringRender.cpp
              ../../Algorithms/library/polygonRender.cpp
              ${PROJECT_SOURCE_DIR}/source-code/exp3_fig15.cpp
              )
target_link_libraries(exp3_fig15 ${HIGIS} ${PNG_LIBRARIES} ${Boost_LIBRARIES} ${MAPNIK} ${MPI_LIBRARIES} -lmpich -lpthread)


add_executable(exp3_fig16  
              ../../Algorithms/library/buildIndex.cpp
              ../../Algorithms/library/indexNode.cpp 
              ../../Algorithms/library/topoJudge.cpp 
              ../../Algorithms/library/pointRender.cpp
              ../../Algorithms/library/linestringRender.cpp
              ../../Algorithms/library/polygonRender.cpp
              ${PROJECT_SOURCE_DIR}/source-code/exp3_fig16.cpp
              )
target_link_libraries(exp3_fig16 ${HIGIS} ${PNG_LIBRARIES} ${Boost_LIBRARIES} ${MAPNIK} ${MPI_LIBRARIES} -lmpich -lpthread)