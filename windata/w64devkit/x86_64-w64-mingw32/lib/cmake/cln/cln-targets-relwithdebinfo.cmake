#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "cln::cln" for configuration "RelWithDebInfo"
set_property(TARGET cln::cln APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(cln::cln PROPERTIES
  IMPORTED_IMPLIB_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libcln.dll.a"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/bin/libcln.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS cln::cln )
list(APPEND _IMPORT_CHECK_FILES_FOR_cln::cln "${_IMPORT_PREFIX}/lib/libcln.dll.a" "${_IMPORT_PREFIX}/bin/libcln.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
