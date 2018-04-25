cmake_minimum_required(VERSION 2.8.7)
if(MSVC)
  # CMake 3.4 introduced a WINDOWS_EXPORT_ALL_SYMBOLS target property that makes it possible to
  # build shared libraries without using the usual declspec() decoration.
  # See: https://blog.kitware.com/create-dlls-on-windows-without-declspec-using-new-cmake-export-all-feature/
  # and https://cmake.org/cmake/help/v3.5/prop_tgt/WINDOWS_EXPORT_ALL_SYMBOLS.html
  # for details.
  cmake_minimum_required(VERSION 3.4)
endif()
if(POLICY CMP0046)
  cmake_policy(SET CMP0046 NEW)
endif()
if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()

project (RawScaler C CXX)
# The version number.
set (RawScaler_VERSION_MAJOR 1)
set (RawScaler_VERSION_MINOR 0)
option(OpenCV_Enable "for hardware design" ON)
option(openmp "for hardware design" ON)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)


if("${CMAKE_GENERATOR}" MATCHES "(Win64|IA64)")
set(arch_hint "x64")
elseif("${CMAKE_GENERATOR_PLATFORM}" MATCHES "ARM64")
set(arch_hint "ARM64")
elseif("${CMAKE_GENERATOR}" MATCHES "ARM")
set(arch_hint "ARM")
elseif("${CMAKE_SIZEOF_VOID_P}" STREQUAL "8")
set(arch_hint "x64")
elseif("$ENV{LIB}" MATCHES "(amd64|ia64)")
set(arch_hint "x64")
endif()

if(NOT arch_hint)
set(arch_hint "x86")
endif()
#
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)
set(RawScaler_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(RawScaler_SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(RawScaler_BIN_DIR ${PROJECT_SOURCE_DIR}/bin)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/RawScaler.config.h.in"
  "${PROJECT_BINARY_DIR}/RawScaler.config.h"
  )
 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")
 
# add the executable
file(GLOB RawScaler_TOOL
    "tools/RawScaler/*.h"
    "tools/RawScaler/*.cpp"
)
add_executable(RawScaler ${RawScaler_TOOL})

file(GLOB RawScaler_SRC
    "include/*.h"
	"include/*.hpp"
    "src/*.cpp"
)

add_library(RawScalerLib ${RawScaler_SRC})
target_link_libraries (RawScaler RawScalerLib)



include_directories(RawScaler ${PROJECT_SOURCE_DIR}/include)
include_directories(RawScalerLib ${PROJECT_SOURCE_DIR}/include)




IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Debug CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)



IF(fixpoint_cut MATCHES ON) 
	add_definitions(-D_enable_fix_point_cut)
ENDIF(fixpoint_cut MATCHES ON) 



IF(OpenCV_Enable MATCHES ON) 
	add_definitions(-DOpenCV)
	add_definitions(-DImagePath="${PROJECT_SOURCE_DIR}/data/kodim19.png")
	
	
ELSE()
	add_definitions(-DImagePath="${PROJECT_SOURCE_DIR}/data/kodim19.bmp")
	
ENDIF() 
add_definitions(-DImagePath2="${PROJECT_SOURCE_DIR}/data/IMG_0040.bmp")
add_definitions(-DImagePath3="${PROJECT_SOURCE_DIR}/data/IMG_0086.bmp")
add_definitions(-DImagePath4="${PROJECT_SOURCE_DIR}/data/kodim01.bmp")

add_definitions(-DOutputDir="${PROJECT_SOURCE_DIR}/out/")

if("${CMAKE_GENERATOR}" MATCHES "Visual Studio 14 2015 Win64")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
	IF(openmp MATCHES ON) 
			add_definitions(-Dopenmp)
	ENDIF(openmp MATCHES ON) 
	message("\n============== Windows " ${arch_hint} " build ===================\n")
	message(STATUS "CMAKE_GENERATOR: ${CMAKE_GENERATOR}")
	set_property( DIRECTORY PROPERTY VS_STARTUP_PROJECT "RawScaler" )
	IF(CMAKE_BUILD_TYPE MATCHES Debug)
		set(OpenCV_DIR OpenCV_DEBUG)
	ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)
	
	IF(CMAKE_BUILD_TYPE MATCHES Release)
		set(OpenCV_DIR OpenCV_RELEASE)
	ENDIF(CMAKE_BUILD_TYPE MATCHES Release)
	
	set(OpenCV_FOUND 1 )

	find_package(OpenCV REQUIRED COMPONENTS core imgproc highgui)
	# OPENCV config
	set(OPENCV_DIR ${CMAKE_CURRENT_LIST_DIR}  CACHE PATH "")
	set(opencv_DIR ${CMAKE_CURRENT_LIST_DIR}  CACHE PATH "")
	set(OpenCV_DIR ${CMAKE_CURRENT_LIST_DIR}  CACHE PATH "")
	set(OpenCV_STATIC OFF CACHE BOOL "")
	# Add OpenCV to search directories
	get_filename_component(_dir ${OpenCV_LIB_PATH} DIRECTORY)
	list(APPEND _directories ${_dir}/bin)
	
	include_directories(${OpenCV_INCLUDE_DIRS})
	target_link_libraries(RawScaler ${OpenCV_LIBS})
	#include_directories(RawScaler "C:/Source Project 2017/RawScaler/OpenCV/x64/vc14/bin")
	#file(GLOB OpenCV_DLLS
	#  "C:/Source Project 2017/RawScaler/OpenCV/x64/vc14/bin/*.dll"
	#)
	#file(COPY ${OpenCV_DLLS} DESTINATION ${RawScaler_BIN_DIR})
	IF(CMAKE_BUILD_TYPE MATCHES Debug)
		
		add_custom_command(TARGET RawScaler PRE_BUILD        # Adds a post-build event to MyTest
		COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
			"${_OpenCV_LIB_PATH}/opencv_highgui2413d.dll"      # <--this is in-file
			$<TARGET_FILE_DIR:RawScaler>)                 # <--this is out-file path
		#target_link_libraries (RawScaler util)
		add_custom_command(TARGET RawScaler PRE_BUILD        # Adds a post-build event to MyTest
		COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
			"${_OpenCV_LIB_PATH}/opencv_imgproc2413d.dll"      # <--this is in-file
			$<TARGET_FILE_DIR:RawScaler>)                 # <--this is out-file path
		add_custom_command(TARGET RawScaler PRE_BUILD        # Adds a post-build event to MyTest
		COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
			"${_OpenCV_LIB_PATH}/opencv_core2413d.dll"      # <--this is in-file
			$<TARGET_FILE_DIR:RawScaler>)                 # <--this is out-file path
	ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)
	
	IF(CMAKE_BUILD_TYPE MATCHES Release)
		add_custom_command(TARGET RawScaler PRE_BUILD        # Adds a post-build event to MyTest
		COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
			"${_OpenCV_LIB_PATH}/opencv_highgui2413.dll"      # <--this is in-file
			$<TARGET_FILE_DIR:RawScaler>)                 # <--this is out-file path
		#target_link_libraries (RawScaler util)
		add_custom_command(TARGET RawScaler PRE_BUILD        # Adds a post-build event to MyTest
		COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
			"${_OpenCV_LIB_PATH}/opencv_imgproc2413.dll"      # <--this is in-file
			$<TARGET_FILE_DIR:RawScaler>)                 # <--this is out-file path
		add_custom_command(TARGET RawScaler PRE_BUILD        # Adds a post-build event to MyTest
		COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
			"${_OpenCV_LIB_PATH}/opencv_core2413.dll"      # <--this is in-file
			$<TARGET_FILE_DIR:RawScaler>)                 # <--this is out-file path
	ENDIF(CMAKE_BUILD_TYPE MATCHES Release)

else()	
	message("\n============== Linux " ${arch_hint} " build ===================\n")
	message(STATUS "CMAKE_GENERATOR: ${CMAKE_GENERATOR}")
	find_package( OpenCV REQUIRED )
	target_link_libraries( RawScaler ${OpenCV_LIBS} )
	#message(FATAL_ERROR "CMAKE_GENERATOR do net match Visual Studio 14 2015 Win64")
	IF(openmp MATCHES ON) 
		find_package(OpenMP)
		if (OPENMP_FOUND)
			set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
			set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
			set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
			#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
			add_definitions(-Dopenmp)
		endif()
	
	ENDIF(openmp MATCHES ON) 

endif()



