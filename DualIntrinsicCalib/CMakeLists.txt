cmake_minimum_required(VERSION 2.8)                         
project( DualIntrinsicCalib )                                                               
find_package( OpenCV REQUIRED )                           
file(GLOB SRC
    "*.h"
    "*.cpp"
) 
add_executable( DualIntrinsicCalib ${SRC})               
target_link_libraries( DualIntrinsicCalib ${OpenCV_LIBS} )  
include_directories($ENV{AC_DIR})

      