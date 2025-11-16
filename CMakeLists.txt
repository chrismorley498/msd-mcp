cmake_minimum_required(VERSION 3.5)

project(msd-mpc)

set(CMAKE_BUILD_TYPE Release)

add_executable(${PROJECT_NAME} main.cpp MPC.hpp)

#matplotlib includes/libs
set(MATPLOTLIB_INCLUDE_DIRS /usr/include/python3.12 /usr/include/python3.12 /home/chris-morley/sandbox/plotting_demo/venv/lib/python3.12/site-packages/numpy/_core/include)
set(MATPLOTLIB_LIBS python3.12)

target_include_directories(${PROJECT_NAME} PUBLIC ${MATPLOTLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${MATPLOTLIB_LIBS} nlopt)
