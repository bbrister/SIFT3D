# FindNIFTI.cmake
#
# Finds the MinGW runtime dependencies:
#
#   -libgcc
#   -libstdc++
#   -libgfortran
#   -libquadmath
#
# The following variables will be exported:
#   MINGW_LIBRARIES - The dependencies

# Standard MinGW locations
set (MINGW_HINTS  "C:/TDM-GCC/")
set (MINGW_SUFFIXES "bin")


# Get the names of the libraries, depending on 32- or 64-bit
if ("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
	message(STATUS "Using 64-bit MinGW...")
	set (MINGW_GCC_LIBRARY_NAME "libgcc_s_seh_64-1")
	set (MINGW_CXX_LIBRARY_NAME "libstdc++_64-6")
	set (MINGW_GFORTRAN_LIBRARY_NAME "libgfortran_64-3")
	set (MINGW_QUADMATH_LIBRARY_NAME "libquadmath_64-0")
elseif ("${CMAKE_SIZEOF_VOID_P}" EQUAL "4")
	message(STATUS "Using 32-bit MinGW...")
	set (MINGW_GCC_LIBRARY_NAME "libgcc_s_dw2-1")
	set (MINGW_CXX_LIBRARY_NAME "libstdc++-6")
	set (MINGW_GFORTRAN_LIBRARY_NAME "libgfortran-3")
	set (MINGW_QUADMATH_LIBRARY_NAME "libquadmath-0")
else ()
        message(FATAL_ERROR "Unrecognized byte width: ${CMAKE_SIZE_OF_VOID_P}")
endif ()


find_library (
	MINGW_GCC_LIBRARY
	NAMES ${MINGW_GCC_LIBRARY_NAME} 
	HINTS ${MINGW_HINTS}
	PATH_SUFFIXES ${MINGW_SUFFIXES} 
	DOC "MinGW C standard library"
)

find_library (
	MINGW_CXX_LIBRARY
	names ${MINGW_CXX_LIBRARY_NAME}
	HINTS ${MINGW_HINTS}
	PATH_SUFFIXES ${MINGW_SUFFIXES} 
	DOC "MinGW C++ standard library"	
)

find_library (
	MINGW_GFORTRAN_LIBRARY
	NAMES ${MINGW_GFORTRAN_LIBRARY_NAME}
	HINTS ${MINGW_HINTS}
	PATH_SUFFIXES ${MINGW_SUFFIXES}
	DOC "MinGW Fortran standard library"
)

find_library (
	MINGW_QUADMATH_LIBRARY
	NAMES ${MINGW_QUADMATH_LIBRARY_NAME}
	HINTS ${MINGW_HINTS}
	PATH_SUFFIXES ${MINGW_SUFFIXES}
	DOC "MinGW quadmath library"
)

set (MINGW_LIBRARIES ${MINGW_GCC_LIBRARY} ${MINGW_CXX_LIBRARY} ${MINGW_GFORTRAN_LIBRARY} ${MINGW_QUADMATH_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set PNG_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MINGW
	"Cannot find package MinGW."
	MINGW_GCC_LIBRARY
	MINGW_CXX_LIBRARY
	MINGW_GFORTRAN_LIBRARY
	MINGW_QUADMATH_LIBRARY
)

# these variables are only visible in 'advanced mode'
MARK_AS_ADVANCED(
	MINGW_GCC_LIBRARY
	MINGW_CXX_LIBRARY
	MINGW_GFORTRAN_LIBRARY
	MINGW_QUADMATH_LIBRARY
)
