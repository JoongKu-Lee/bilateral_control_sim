#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <iostream>
#include <math.h>
#include <Eigen/Dense>

class DynamicsSolver
{
private:
    double L01_z = 0.07;
    double L12_x = 0.08;
    double L12_z = 0.12;
    double L23_y = 0.43;
    double L34_y = 0.335;
    double L34_z = -0.08;
    double L1_xcom = 0.008355;
    double L1_zcom = 0.102659;
    double L2_ycom = 0.213502;
    double L2_zcom = 0.0919281;
    double L3_ycom = 0.120975;
    double L3_zcom = -0.0755858;
    double m1 = 2.75458;
    double m2 = 8.143;
    double m3 = 5.14886;
    double I1_xx = 0.0109697;
    double I1_yy = 0.011856;
    double I1_zz = 0.00604953;
    double I1_xy = -2.08287e-08;
    double I1_yz = 2.03974e-09;
    double I1_zx = -0.000398309;
    double I2_xx = 0.25374;
    double I2_yy = 0.0212824;
    double I2_zz = 0.247523;
    double I2_xy = 1.41285e-07;
    double I2_yz = 0.00499483;
    double I2_zx = -7.13886e-08;
    double I3_xx = 0.075434;
    double I3_yy = 0.010425;
    double I3_zz = 0.0744749;
    double I3_xy = 2.67415e-08;
    double I3_yz = 0.00274946;
    double I3_zx = 2.13561e-09;
    double g = 9.81;
    
public:
    Eigen::Matrix3d Mass(const Eigen::Vector3d &q);
    Eigen::Matrix3d Coriolis(const Eigen::Vector3d &q, const Eigen::Vector3d &dq);
    Eigen::Vector3d Gravity(const Eigen::Vector3d &q);
};



#endif