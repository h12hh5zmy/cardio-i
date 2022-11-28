cmake_minimum_required(VERSION 3.2)
project(ExtractHolterImagesMFER CXX)

find_package( OpenCV  REQUIRED )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++14")

link_libraries(pthread gsl gslcblas)

include_directories(include HolterDataManager qrs_detector/cpp)
add_executable(ExtractHolterImages ExtractImagesPolylines.cpp HolterDataManager/HolterDataManager.cpp HolterDataManager/HolterData.cpp qrs_detector/cpp/QRSdetector.cpp)

target_link_libraries(ExtractHolterImages ${OpenCV_LIBS})
