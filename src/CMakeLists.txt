set(mc_zmp_area_src
  mc_zmp_area.cpp
)

add_library(mc_zmp_area SHARED ${mc_zmp_area_src})

target_link_libraries(mc_zmp_area mc_rtc::mc_control)
target_link_libraries(mc_zmp_area eigen-lssol::eigen-lssol)

pkg_config_use_dependency(mc_zmp_area nlopt)

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION include/mc_zmp_area
    FILES_MATCHING PATTERN "*.h*")

install(TARGETS mc_zmp_area
  EXPORT "${TARGETS_EXPORT_NAME}"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )



