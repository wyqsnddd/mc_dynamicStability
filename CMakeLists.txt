cmake_minimum_required(VERSION 3.8)

include(cmake/base.cmake)

set(PROJECT_NAME mc_impact_predictor)

set(PROJECT_DESCRIPTION "Calculate the multi-contact ZMP area.")
set(PROJECT_URL https://gite.lirmm.fr/yuquan/mc_zmp_area)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CXX_DISABLE_WERROR 1)
set(CMAKE_CXX_STANDARD 11)
setup_project()

find_package(mc_rtc REQUIRED)

add_required_dependency(nlopt)

include_directories(include)
add_subdirectory(src)

setup_project_finalize()
setup_project_package_finalize()

