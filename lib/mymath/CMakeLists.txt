set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(UNIX)
  find_library(ginac ginac REQUIRED)
  find_library(cln cln REQUIRED)
endif()

add_library(mymath STATIC mymath.cpp)
target_include_directories(mymath PUBLIC include)
target_link_libraries(mymath cln ginac)
