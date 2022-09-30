cmake_minimum_required(VERSION 3.21)
project(eq36 Fortran)

set(${PROJECT_NAME}_MAJOR 80)
set(${PROJECT_NAME}_MINOR 1)
set(${PROJECT_NAME}_PATCH 0)
set(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_MAJOR}.${${PROJECT_NAME}_MINOR}.${${PROJECT_NAME}_PATCH})

message(STATUS "eq36 version: ${${PROJECT_NAME}_VERSION}")

set(CMAKE_FORTRAN_FLAGS_DEBUG "${CMAKE_FORTRAN_FLAGS_DEBUG} -g -pg")

include_directories(include include/eqlib include/xcon3 include/xcon6)
add_subdirectory(src/eqlib)
add_subdirectory(src/eqlibu)
add_subdirectory(src/eqlibg)
add_subdirectory(src/eqpt)
add_subdirectory(src/eq3nr)
add_subdirectory(src/eq6)
add_subdirectory(src/xcon3)
add_subdirectory(src/xcon6)

install(DIRECTORY include/ DESTINATION include)
