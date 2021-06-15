#include <iostream>
#include <Eigen/Dense>
#include "myrobot_control/Dynamics.h"

using namespace std;

int main()
{
    Eigen::Vector3d asd(-0.1, 0.2, 0.3);

    Dynamics robot;
    Eigen::Matrix3d Mass = robot.Mass(asd);
    Eigen::Matrix3d Coriolis = robot.Coriolis(asd, asd);
    Eigen::Vector3d Gravity = robot.Gravity(asd);

    cout<<Mass<<endl;
    cout<<Coriolis<<endl;

    cout<<Mass - 2 * Coriolis <<endl;
    
    // cout<<Gravity.transpose()<<endl;


    return 0;
}