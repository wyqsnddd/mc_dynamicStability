set(mc_zmp_area_src
  McContact.cpp
  McZMPArea.cpp
  Utils.cpp
)

add_library(${PROJECT_NAME} SHARED ${mc_zmp_area_src})

target_link_libraries(${PROJECT_NAME} mc_rtc::mc_control)
target_link_libraries(${PROJECT_NAME} eigen-lssol::eigen-lssol)
target_link_libraries(${PROJECT_NAME} Stabiliplus::StabilityPolytope)
target_link_libraries(${PROJECT_NAME} Stabiliplus::ProblemDescriptor)
target_link_libraries(${PROJECT_NAME} Stabiliplus::wrappers)
target_link_libraries(${PROJECT_NAME} qhullcpp qhullstatic_r)

pkg_config_use_dependency(${PROJECT_NAME} nlopt)

include(GNUInstallDirs)

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h*")

install(TARGETS ${PROJECT_NAME} 
  EXPORT "${TARGETS_EXPORT_NAME}"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )



