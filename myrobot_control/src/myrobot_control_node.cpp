#include <iostream>
#include <Eigen/Dense>
#include "myrobot_control/DynamicsSolver.h"
#include <ros/ros.h>
#include <ros/publisher.h>
#include <ros/subscriber.h>
#include <sensor_msgs/JointState.h>
#include <std_msgs/Float64.h>

using namespace std;
using namespace Eigen;

ros::Publisher myrobot1_joint1_effort_pub;
ros::Publisher myrobot1_joint2_effort_pub;
ros::Publisher myrobot1_joint3_effort_pub;
ros::Publisher myrobot2_joint1_effort_pub;
ros::Publisher myrobot2_joint2_effort_pub;
ros::Publisher myrobot2_joint3_effort_pub;


void callback_myrobot_joint_state(const sensor_msgs::JointStateConstPtr &msg)
{
    Vector3d q(msg->position[0], msg->position[1], msg->position[2]);
    Vector3d dq(msg->velocity[0], msg->velocity[1], msg->velocity[2]);

    DynamicsSolver DS;
    
    Matrix3d Mass = DS.Mass(q);
    Matrix3d Coriolis = DS.Coriolis(q, dq);
    Vector3d Gravity = DS.Gravity(q);

    cout<<Gravity.transpose()<<endl;
    
    std_msgs::Float64 joint1_effort;
    std_msgs::Float64 joint2_effort;
    std_msgs::Float64 joint3_effort;

    joint1_effort.data = Gravity(0);
    joint2_effort.data = Gravity(1);
    joint3_effort.data = Gravity(2);

    myrobot1_joint1_effort_pub.publish(joint1_effort);
    myrobot1_joint2_effort_pub.publish(joint2_effort);
    myrobot1_joint3_effort_pub.publish(joint3_effort);
}

void callback_myrobot2_joint_state(const sensor_msgs::JointStateConstPtr &msg)
{
    Vector3d q(msg->position[0], msg->position[1], msg->position[2]);
    Vector3d dq(msg->velocity[0], msg->velocity[1], msg->velocity[2]);

    DynamicsSolver DS;
    
    Matrix3d Mass = DS.Mass(q);
    Matrix3d Coriolis = DS.Coriolis(q, dq);
    Vector3d Gravity = DS.Gravity(q);

    // cout<<Gravity.transpose()<<endl;
    
    std_msgs::Float64 joint1_effort;
    std_msgs::Float64 joint2_effort;
    std_msgs::Float64 joint3_effort;

    joint1_effort.data = Gravity(0);
    joint2_effort.data = Gravity(1);
    joint3_effort.data = Gravity(2);

    // myrobot2_joint1_effort_pub.publish(joint1_effort);
    // myrobot2_joint2_effort_pub.publish(joint2_effort);
    // myrobot2_joint3_effort_pub.publish(joint3_effort);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "myrobot_control_node");

    ros::NodeHandle nh;
    ros::NodeHandle pnh("~");

    ros::Subscriber myrobot1_joint_state_sub;
    myrobot1_joint_state_sub = nh.subscribe("/myrobot/joint_states", 3, callback_myrobot_joint_state);

    myrobot1_joint1_effort_pub = nh.advertise<std_msgs::Float64>("/myrobot/joint1_effort_controller/command", 3);
    myrobot1_joint2_effort_pub = nh.advertise<std_msgs::Float64>("/myrobot/joint2_effort_controller/command", 3);
    myrobot1_joint3_effort_pub = nh.advertise<std_msgs::Float64>("/myrobot/joint3_effort_controller/command", 3);

    ros::Subscriber myrobot2_joint_state_sub;
    myrobot2_joint_state_sub = nh.subscribe("/myrobot2/joint_states", 3, callback_myrobot2_joint_state);

    myrobot2_joint1_effort_pub = nh.advertise<std_msgs::Float64>("/myrobot2/joint1_effort_controller/command", 3);
    myrobot2_joint2_effort_pub = nh.advertise<std_msgs::Float64>("/myrobot2/joint2_effort_controller/command", 3);
    myrobot2_joint3_effort_pub = nh.advertise<std_msgs::Float64>("/myrobot2/joint3_effort_controller/command", 3);

    ros::spin();

    return 0;
}