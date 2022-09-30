cmake_minimum_required(VERSION 3.21)
project(xcon3 Fortran)

set(${PROJECT_NAME}_MAJOR 80)
set(${PROJECT_NAME}_MINOR 1)
set(${PROJECT_NAME}_PATCH 0)
set(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_MAJOR}.${${PROJECT_NAME}_MINOR}.${${PROJECT_NAME}_PATCH})

message(STATUS "${PROJECT_NAME} version: ${${PROJECT_NAME}_VERSION}")

set(CMAKE_FORTRAN_FLAGS_DEBUG "${CMAKE_FORTRAN_FLAGS_DEBUG} -g -pg")

add_subdirectory(src)
add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})
target_link_libraries(${PROJECT_NAME} eqlibu_static)

install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME} DESTINATION bin)
