cmake_minimum_required(VERSION 3.21)
project(eq6 Fortran)

set(${PROJECT_NAME}_MAJOR 80)
set(${PROJECT_NAME}_MINOR 1)
set(${PROJECT_NAME}_PATCH 0)
set(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_MAJOR}.${${PROJECT_NAME}_MINOR}.${${PROJECT_NAME}_PATCH})

message(STATUS "${PROJECT_NAME} version: ${${PROJECT_NAME}_VERSION}")

set(CMAKE_FORTRAN_FLAGS_DEBUG "${CMAKE_FORTRAN_FLAGS_DEBUG} -g -pg")

add_subdirectory(src)

add_library(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME} VERSION ${${PROJECT_NAME}_VERSION})
target_link_libraries(${PROJECT_NAME} eqlib eqlibu eqlibg)

add_library(${PROJECT_NAME}_static ${${PROJECT_NAME}_SOURCES})
set_target_properties(${PROJECT_NAME}_static PROPERTIES OUTPUT_NAME ${PROJECT_NAME} VERSION ${${PROJECT_NAME}_VERSION})

add_executable(${PROJECT_NAME}_exe src/eq6.f)
set_target_properties(${PROJECT_NAME}_exe PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}_exe ${PROJECT_NAME}_static eqlib_static eqlibu_static eqlibg_static)

install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME} DESTINATION lib)
install(TARGETS ${PROJECT_NAME}_static EXPORT ${PROJECT_NAME} DESTINATION lib)
install(TARGETS ${PROJECT_NAME}_exe EXPORT ${PROJECT_NAME} DESTINATION bin)
