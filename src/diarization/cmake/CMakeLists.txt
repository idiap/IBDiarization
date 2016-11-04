project(idiap-aib-diarization)
cmake_minimum_required(VERSION 2.8)

set(AIBROOT ../../)
set(CMAKE_CXX_FLAGS " -Wall -fopenmp ")

file(
    GLOB_RECURSE
	source_files
	${AIBROOT}/diarization/src/*
	${AIBROOT}/aIB/src/aIB.cpp
	${AIBROOT}/aIB/src/aibfunctions.cpp
	${AIBROOT}/aIB/src/fileIO.cpp
	${AIBROOT}/aIB/include/*
	${AIBROOT}/utils/src/*
	${AIBROOT}/utils/include/*
	${AIBROOT}/aibfeat/src/*
	${AIBROOT}/aibfeat/include/*
	${AIBROOT}/realign/src/*
	${AIBROOT}/realign/include/*
)


find_package(Boost COMPONENTS system filesystem program_options REQUIRED)

include_directories(${Boost_INCLUDE_DIRS} : ${AIBROOT}/aIB/include/ : ${AIBROOT}/utils/include/ : ${AIBROOT}/aibfeat/include/ : ${AIBROOT}/realign/include/ )

add_executable(
  diarizeme
  ${source_files}
)

target_link_libraries( diarizeme ${Boost_LIBRARIES} )
install (TARGETS diarizeme DESTINATION bin)
