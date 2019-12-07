cmake_minimum_required(VERSION 3.8)

project(mc_dynamicStability)

set(PROJECT_URL https://gite.lirmm.fr/jrl-collaborators/mc_zmp_area)
set(PROJECT_NAME mc_dynamicStability)
set(PROJECT_DESCRIPTION "Calculate the multi-contact ZMP area.")


include(cmake/base.cmake)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CXX_DISABLE_WERROR 1)
set(CMAKE_CXX_STANDARD 11)
setup_project()

find_package(mc_rtc REQUIRED)
find_package(Stabiliplus REQUIRED)

add_required_dependency(nlopt)

include_directories(include)
add_subdirectory(src)

setup_project_finalize()
setup_project_package_finalize()

