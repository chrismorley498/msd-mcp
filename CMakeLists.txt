cmake_minimum_required(VERSION 3.5)

project(msd-mpc)

# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

add_library(${PROJECT_NAME} include/MPC.hpp)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX) 

#matplotlib includes/libs
set(MATPLOTLIB_INCLUDE_DIRS /usr/include/python3.12 /usr/include/python3.12 /home/chris-morley/sandbox/plotting_demo/venv/lib/python3.12/site-packages/numpy/_core/include)
set(MATPLOTLIB_LIBS python3.12)

target_include_directories(${PROJECT_NAME} PUBLIC ${MATPLOTLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${MATPLOTLIB_LIBS} nlopt)

#MPC test
add_executable(mpc_test test/mpc_test.cpp)
target_include_directories(mpc_test PUBLIC include)
target_link_libraries(mpc_test ${PROJECT_NAME})

#Simulation test
add_executable(simulation_test test/simulation_test.cpp)
target_include_directories(simulation_test PUBLIC include)
target_link_libraries(simulation_test ${PROJECT_NAME})