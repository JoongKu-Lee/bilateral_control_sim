cmake_minimum_required(VERSION 3.0.2)
project(myrobot_control)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
)

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES virtual_master_device
  CATKIN_DEPENDS roscpp rospy std_msgs 
#  DEPENDS system_lib
)

find_package(Eigen3 REQUIRED)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}_lib SHARED
        include/myrobot_control/Dynamics.cpp
)

add_executable(test_node src/test.cpp)
target_link_libraries(test_node
  ${catkin_LIBRARIES}
  ${PROJECT_NAME}_lib
)