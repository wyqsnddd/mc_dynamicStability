set(mc_zmp_area_src
  McContact.cpp
  McZMPArea.cpp
  McComArea.cpp
  Utils.cpp
)

add_library(${PROJECT_NAME} SHARED ${mc_zmp_area_src})

target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} PUBLIC mc_rtc::mc_control)
target_link_libraries(${PROJECT_NAME} PUBLIC eigen-lssol::eigen-lssol)
target_link_libraries(${PROJECT_NAME} PUBLIC Stabiliplus::StabilityPolytope)
target_link_libraries(${PROJECT_NAME} PUBLIC Stabiliplus::ProblemDescriptor)
target_link_libraries(${PROJECT_NAME} PUBLIC Stabiliplus::wrappers)
target_link_libraries(${PROJECT_NAME} PRIVATE Qhull::qhullcpp)
target_link_libraries(${PROJECT_NAME} PUBLIC nlopt)

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



