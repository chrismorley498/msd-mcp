cmake_minimum_required(VERSION 3.5)

project(msd-mpc)

set(pinocchio_DIR "/opt/openrobots/lib/cmake/pinocchio")
set(hpp-fcl_DIR "/opt/openrobots/lib/cmake/hpp-fcl")
find_package(pinocchio REQUIRED)
find_package(hpp-fcl REQUIRED)

# Force Miniconda Boost
set(BOOST_ROOT /home/chris-morley/miniconda3)
set(BOOST_INCLUDEDIR /home/chris-morley/miniconda3/include)
set(BOOST_LIBRARYDIR /home/chris-morley/miniconda3/lib)
find_package(Boost 1.88 REQUIRED COMPONENTS system filesystem) # add your needed components
include_directories(BEFORE ${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

message(STATUS "Boost version: ${Boost_VERSION}")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth=1024")



# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

add_library(${PROJECT_NAME} include/MPC.hpp)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX) 

#matplotlib includes/libs
set(MATPLOTLIB_INCLUDE_DIRS /usr/include/python3.12 /usr/include/python3.12 /home/chris-morley/sandbox/plotting_demo/venv/lib/python3.12/site-packages/numpy/_core/include)
set(MATPLOTLIB_LIBS python3.12)

target_include_directories(${PROJECT_NAME} PUBLIC ${MATPLOTLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${MATPLOTLIB_LIBS} nlopt pinocchio::pinocchio )

#MPC test
add_executable(mpc_test test/mpc_test.cpp)
target_include_directories(mpc_test PUBLIC include)
target_link_libraries(mpc_test ${PROJECT_NAME})

#Simulation test
add_executable(simulation_test test/simulation_test.cpp)
target_include_directories(simulation_test PUBLIC include)
target_link_libraries(simulation_test ${PROJECT_NAME})


#Three DoF test
add_executable(three_dof_test test/three_dof_sim_test.cpp)
target_include_directories(three_dof_test PUBLIC include)
target_link_libraries(three_dof_test ${PROJECT_NAME})