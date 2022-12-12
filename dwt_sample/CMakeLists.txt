cmake_minimum_required(VERSION 3.2)
project(dwt_sample CXX)

include_directories(include HolterDataManager HolterDataManager/include qrs_detector/cpp)

add_executable(${PROJECT_NAME} dwt_sample.cpp HolterDataManager/HolterDataManager.cpp HolterDataManager/HolterData.cpp qrs_detector/cpp/QRSdetector.cpp)

# wavelet
include_directories("C:/local/dwt_sample/wavelet/include")
target_link_libraries(${PROJECT_NAME} "C:/local/dwt_sample/wavelet/wauxlib.lib")
target_link_libraries(${PROJECT_NAME} "C:/local/dwt_sample/wavelet/wavelib.lib")

target_link_libraries(${PROJECT_NAME})

