cmake_minimum_required(VERSION 3.16.0)

project(Mesh2Sdf)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(VTK REQUIRED)

include(${VTK_USE_FILE})

add_executable(${PROJECT_NAME} Src/main.cpp Src/AabbBox.h Src/AabbBox.cpp Src/SignedDistanceField.h Src/SignedDistanceField.cpp)

target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})