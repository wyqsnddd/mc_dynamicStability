set(mc_zmp_area_src
  contact.cpp
  mc_zmp_area.cpp
)

add_library(${PROJECT_NAME} SHARED ${mc_zmp_area_src})

target_link_libraries(${PROJECT_NAME} mc_rtc::mc_control)
target_link_libraries(${PROJECT_NAME} eigen-lssol::eigen-lssol)
target_link_libraries(${PROJECT_NAME} Stabiliplus::StabilityPolytope)

pkg_config_use_dependency(${PROJECT_NAME} nlopt)

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION include/${PROJECT_NAME}
    FILES_MATCHING PATTERN "*.h*")

install(TARGETS ${PROJECT_NAME} 
  EXPORT "${TARGETS_EXPORT_NAME}"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )



