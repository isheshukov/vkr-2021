include(CMakeFindDependencyMacro)

set(CLN_USE_GMP ON)

if (CLN_USE_GMP)
	list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
	find_package(GMP REQUIRED MODULE)
	list(REMOVE_AT CMAKE_MODULE_PATH -1)
endif()

if (NOT TARGET cln::cln)
	include("${CMAKE_CURRENT_LIST_DIR}/cln-targets.cmake")
endif()

set(cln_LIBRARIES cln::cln)
