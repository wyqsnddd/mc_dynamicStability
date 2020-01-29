cmake_minimum_required(VERSION 3.1)

set(PROJECT_NAME mc_dynamicStability)
set(PROJECT_VERSION 0.1.0)
set(PROJECT_DESCRIPTION "Calculate the multi-contact ZMP area.")
set(PROJECT_URL https://gite.lirmm.fr/mc-controllers/mc_dynamicStability)
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CXX_DISABLE_WERROR 1)
set(CMAKE_CXX_STANDARD 11)
set(PROJECT_DEBUG_POSTFIX "_d")
set(DOXYGEN_USE_MATHJAX "YES")

include(cmake/base.cmake)

project(${PROJECT_NAME} LANGUAGES CXX VERSION ${PROJECT_VERSION})

add_project_dependency(mc_rtc REQUIRED)
add_project_dependency(eigen-lssol REQUIRED)
add_project_dependency(Stabiliplus REQUIRED)
add_project_dependency(Qhull 7.3.2 REQUIRED)

add_subdirectory(src)
