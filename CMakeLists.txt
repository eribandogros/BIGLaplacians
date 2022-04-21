cmake_minimum_required(VERSION 3.1)

project(EulerianMesh3D)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

find_package(Matlab)
find_package(Threads)

find_library(LIB_MATLAB_DATA_ARRAY 
	NAMES MatlabDataArray
	PATHS ${Matlab_ROOT_DIR}/extern/bin/maci64
	NO_DEFAULT_PATH)

find_library(LIB_MATLAB_ENGINE
	NAMES MatlabEngine
	PATHS ${Matlab_ROOT_DIR}/extern/bin/maci64
        NO_DEFAULT_PATH)

include_directories(${Matlab_INCLUDE_DIRS})
include_directories(eigen)

file(GLOB Headers "src/*.h")
file(GLOB Sources "src/*.cpp")

add_executable(NormalEulerianMesh3D
	${Headers} 
	${Sources})

target_link_libraries(NormalEulerianMesh3D
	${LIB_MATLAB_DATA_ARRAY}
	${LIB_MATLAB_ENGINE}
	Threads::Threads 	
)	
