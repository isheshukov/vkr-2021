cmake_minimum_required(VERSION 3.5)

project(tropical LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(tropical
  main.cpp
)
target_link_libraries(tropical PUBLIC structopt mymath)

if(WIN32)
add_custom_command(TARGET tropical POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${PROJECT_SOURCE_DIR}/dlls"
        $<TARGET_FILE_DIR:tropical>)
endif()
