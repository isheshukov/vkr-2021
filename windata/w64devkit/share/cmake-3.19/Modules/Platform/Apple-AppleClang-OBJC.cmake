include(Platform/Apple-Clang-OBJC)
if(NOT CMAKE_OBJC_COMPILER_VERSION VERSION_LESS 4.2)
  set(CMAKE_OBJC_SYSTEM_FRAMEWORK_SEARCH_FLAG "-iframework ")
else()
  unset(CMAKE_OBJC_SYSTEM_FRAMEWORK_SEARCH_FLAG)
endif()
