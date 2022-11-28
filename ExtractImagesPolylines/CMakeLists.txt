cmake_minimum_required(VERSION 3.2)
project(ExtractHolterImages CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++14")

include_directories(include HolterDataManager HolterDataManager/include qrs_detector/cpp)

add_executable(${PROJECT_NAME} ExtractImagesPolylines.cpp HolterDataManager/HolterDataManager.cpp HolterDataManager/HolterData.cpp qrs_detector/cpp/QRSdetector.cpp)

# GSL
include_directories("C:/local/ExtractImagesPolylines/GSL/include")
target_link_libraries(${PROJECT_NAME} "C:/local/ExtractImagesPolylines/GSL/lib/gsl/gsl.lib")

# OpenCV
set(OpenCV_DIR "C:/opencv-4.5.4/install")
find_package(OpenCV)

target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})

