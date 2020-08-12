project (sarw)
cmake_minimum_required(VERSION 2.8)

set(JDFTX_SRC "<full-path-to-jdftx-src>" CACHE STRING "Full path to the jdftx source root")
set(JDFTX_BUILD "<full-path-to-jdftx-build>" CACHE STRING "Full path to the jdftx build directory")
set(CMAKE_MODULE_PATH "${JDFTX_SRC}/CMake-Modules/")

#Package configuration:
set(CPACK_PACKAGE_NAME "sarw")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "2")
set(CPACK_PACKAGE_VERSION_PATCH "alpha")
set(VERSION_STRING "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

#Get revision information
find_package(Git)
if(GIT_FOUND)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} rev-parse --is-inside-work-tree
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		OUTPUT_VARIABLE IS_GIT
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if("true" STREQUAL "${IS_GIT}")
		execute_process(
			COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			OUTPUT_VARIABLE GIT_HASH
			OUTPUT_STRIP_TRAILING_WHITESPACE)
		message(STATUS "Git revision hash: ${GIT_HASH}")
	endif()
endif()

#Whether MPI is on or off should match the corresponding JDFTx build:
option(EnableMPI "Use MPI parallelization (in addition to threads / gpu)" ON)
if(EnableMPI)
	find_package(MPI REQUIRED)
	include_directories(${MPI_INCLUDE_PATH})
	add_definitions("-DMPI_ENABLED")
endif()

configure_file(${CMAKE_SOURCE_DIR}/config.in.h ${CMAKE_BINARY_DIR}/config.h)
include_directories("${JDFTX_SRC}")
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR})
add_definitions(-Wall -g -O3 -std=c++0x ${MPI_CXX_COMPILE_FLAGS})

#Profiling (optional):
option(EnableProfiling "Optionally enable profiling to determine ideal functions to optimize (for developers)")
if(EnableProfiling)
	add_definitions("-DENABLE_PROFILING")
endif()

find_library(JDFTX_LIBRARY jdftx PATHS ${JDFTX_BUILD} NO_DEFAULT_PATH)

#Executable:
file(GLOB sarwSources *.cpp)
add_executable(sarw ${sarwSources})
target_link_libraries(sarw ${JDFTX_LIBRARY})
