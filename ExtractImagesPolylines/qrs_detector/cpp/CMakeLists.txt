cmake_minimum_required(VERSION 3.20)
project(qrs_detector)
set(CMAKE_CXX_STANDARD 14)

add_library(qrs_detector STATIC QRSdetector.cpp)
install(TARGETS qrs_detector LIBRARY DESTINATION lib)
install(FILES
	QRSdetector.h
	DESTINATION include)
