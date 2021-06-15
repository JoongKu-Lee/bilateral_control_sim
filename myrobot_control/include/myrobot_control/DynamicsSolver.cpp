#include "DynamicsSolver.h"

Eigen::Matrix3d DynamicsSolver::Mass(const Eigen::Vector3d &q)
{
    Eigen::Matrix3d mass;

    mass(0,0) = I2_xx/2 + I3_xx/2 + I2_yy/2 + I1_zz + I3_yy/2 - (I3_xx*cos(2*q(1) + 2*q(2)))/2 + (I3_yy*cos(2*q(1) + 2*q(2)))/2 + I3_xy*sin(2*q(1) + 2*q(2)) + pow(L12_x, 2)*m2 + pow(L12_x,2)*m3 + (pow(L23_y,2)*m3)/2 + pow(L1_xcom,2)*m1 + (pow(L2_ycom,2)*m2)/2 + pow(L2_zcom,2)*m2 + (pow(L3_ycom,2)*m3)/2 + pow(L3_zcom,2)*m3 - (I2_xx*cos(2*q(1)))/2 + (I2_yy*cos(2*q(1)))/2 + I2_xy*sin(2*q(1)) + 2*L12_x*L2_zcom*m2 + 2*L12_x*L3_zcom*m3 - (pow(L23_y,2)*m3*cos(2*q(1)))/2 - (pow(L2_ycom,2)*m2*cos(2*q(1)))/2 - (pow(L3_ycom,2)*m3*cos(2*q(1) + 2*q(2)))/2 + L23_y*L3_ycom*m3*cos(q(2)) - L23_y*L3_ycom*m3*cos(2*q(1) + q(2));
    mass(0,1) = I3_yz*cos(q(1) + q(2)) + I3_zx*sin(q(1) + q(2)) + I2_yz*cos(q(1)) + I2_zx*sin(q(1)) - L12_x*L3_ycom*m3*cos(q(1) + q(2)) - L3_ycom*L3_zcom*m3*cos(q(1) + q(2)) - L12_x*L23_y*m3*cos(q(1)) - L12_x*L2_ycom*m2*cos(q(1)) - L23_y*L3_zcom*m3*cos(q(1)) - L2_ycom*L2_zcom*m2*cos(q(1));
    mass(0,2) = I3_yz*cos(q(1) + q(2)) + I3_zx*sin(q(1) + q(2)) - L12_x*L3_ycom*m3*cos(q(1) + q(2)) - L3_ycom*L3_zcom*m3*cos(q(1) + q(2));
    mass(1,0) = I3_yz*cos(q(1) + q(2)) + I3_zx*sin(q(1) + q(2)) + I2_yz*cos(q(1)) + I2_zx*sin(q(1)) - L12_x*L3_ycom*m3*cos(q(1) + q(2)) - L3_ycom*L3_zcom*m3*cos(q(1) + q(2)) - L12_x*L23_y*m3*cos(q(1)) - L12_x*L2_ycom*m2*cos(q(1)) - L23_y*L3_zcom*m3*cos(q(1)) - L2_ycom*L2_zcom*m2*cos(q(1));
    mass(1,1) = m3*pow(L23_y,2) + 2*m3*cos(q(2))*L23_y*L3_ycom + m2*pow(L2_ycom,2) + m3*pow(L3_ycom,2) + I2_zz + I3_zz;
    mass(1,2) = m3*pow(L3_ycom,2) + L23_y*m3*cos(q(2))*L3_ycom + I3_zz;
    mass(2,0) = I3_yz*cos(q(1) + q(2)) + I3_zx*sin(q(1) + q(2)) - L12_x*L3_ycom*m3*cos(q(1) + q(2)) - L3_ycom*L3_zcom*m3*cos(q(1) + q(2));
    mass(2,1) = m3*pow(L3_ycom,2) + L23_y*m3*cos(q(2))*L3_ycom + I3_zz;
    mass(2,2) = m3*pow(L3_ycom,2) + I3_zz;

    return mass;
}

Eigen::Matrix3d DynamicsSolver::Coriolis(const Eigen::Vector3d &q, const Eigen::Vector3d &dq)
{
    Eigen::Matrix3d coriolis;
    
    coriolis(0,0) = dq(1)*((m3*sin(2*q(1))*pow(L23_y,2))/2 + m3*sin(2*q(1) + q(2))*L23_y*L3_ycom + (m2*sin(2*q(1))*pow(L2_ycom,2))/2 + (m3*sin(2*q(1) + 2*q(2))*pow(L3_ycom,2))/2 + I3_xy*cos(2*q(1) + 2*q(2)) + (I3_xx*sin(2*q(1) + 2*q(2)))/2 - (I3_yy*sin(2*q(1) + 2*q(2)))/2 + I2_xy*cos(2*q(1)) + (I2_xx*sin(2*q(1)))/2 - (I2_yy*sin(2*q(1)))/2) + dq(2)*(I3_xy*cos(2*q(1) + 2*q(2)) + (I3_xx*sin(2*q(1) + 2*q(2)))/2 - (I3_yy*sin(2*q(1) + 2*q(2)))/2 + (pow(L3_ycom,2)*m3*sin(2*q(1) + 2*q(2)))/2 - (L23_y*L3_ycom*m3*sin(q(2)))/2 + (L23_y*L3_ycom*m3*sin(2*q(1) + q(2)))/2);
    coriolis(0,1) = I2_xy*dq(0)*cos(2*q(1)) - I2_yz*dq(1)*sin(q(1)) + (I2_xx*dq(0)*sin(2*q(1)))/2 - (I2_yy*dq(0)*sin(2*q(1)))/2 + I3_xy*dq(0)*cos(2*q(1) + 2*q(2)) + (I3_xx*dq(0)*sin(2*q(1) + 2*q(2)))/2 - (I3_yy*dq(0)*sin(2*q(1) + 2*q(2)))/2 + I3_zx*dq(1)*cos(q(1) + q(2)) + I3_zx*dq(2)*cos(q(1) + q(2)) - I3_yz*dq(1)*sin(q(1) + q(2)) - I3_yz*dq(2)*sin(q(1) + q(2)) + I2_zx*dq(1)*cos(q(1)) + (pow(L23_y,2)*m3*dq(0)*sin(2*q(1)))/2 + (pow(L2_ycom,2)*m2*dq(0)*sin(2*q(1)))/2 + (pow(L3_ycom,2)*m3*dq(0)*sin(2*q(1) + 2*q(2)))/2 + L12_x*L3_ycom*m3*dq(1)*sin(q(1) + q(2)) + L12_x*L3_ycom*m3*dq(2)*sin(q(1) + q(2)) + L3_ycom*L3_zcom*m3*dq(1)*sin(q(1) + q(2)) + L3_ycom*L3_zcom*m3*dq(2)*sin(q(1) + q(2)) + L12_x*L23_y*m3*dq(1)*sin(q(1)) + L12_x*L2_ycom*m2*dq(1)*sin(q(1)) + L23_y*L3_zcom*m3*dq(1)*sin(q(1)) + L2_ycom*L2_zcom*m2*dq(1)*sin(q(1)) + L23_y*L3_ycom*m3*dq(0)*sin(2*q(1) + q(2));
    coriolis(0,2) = I3_xy*dq(0)*cos(2*q(1) + 2*q(2)) + (I3_xx*dq(0)*sin(2*q(1) + 2*q(2)))/2 - (I3_yy*dq(0)*sin(2*q(1) + 2*q(2)))/2 + I3_zx*dq(1)*cos(q(1) + q(2)) + I3_zx*dq(2)*cos(q(1) + q(2)) - I3_yz*dq(1)*sin(q(1) + q(2)) - I3_yz*dq(2)*sin(q(1) + q(2)) + (pow(L3_ycom,2)*m3*dq(0)*sin(2*q(1) + 2*q(2)))/2 + L12_x*L3_ycom*m3*dq(1)*sin(q(1) + q(2)) + L12_x*L3_ycom*m3*dq(2)*sin(q(1) + q(2)) + L3_ycom*L3_zcom*m3*dq(1)*sin(q(1) + q(2)) + L3_ycom*L3_zcom*m3*dq(2)*sin(q(1) + q(2)) - (L23_y*L3_ycom*m3*dq(0)*sin(q(2)))/2 + (L23_y*L3_ycom*m3*dq(0)*sin(2*q(1) + q(2)))/2;
    coriolis(1,0) = -dq(0)*((m3*sin(2*q(1))*pow(L23_y,2))/2 + m3*sin(2*q(1) + q(2))*L23_y*L3_ycom + (m2*sin(2*q(1))*pow(L2_ycom,2))/2 + (m3*sin(2*q(1) + 2*q(2))*pow(L3_ycom,2))/2 + I3_xy*cos(2*q(1) + 2*q(2)) + (I3_xx*sin(2*q(1) + 2*q(2)))/2 - (I3_yy*sin(2*q(1) + 2*q(2)))/2 + I2_xy*cos(2*q(1)) + (I2_xx*sin(2*q(1)))/2 - (I2_yy*sin(2*q(1)))/2);
    coriolis(1,1) = -L23_y*L3_ycom*m3*dq(2)*sin(q(2));
    coriolis(1,2) = -L23_y*L3_ycom*m3*sin(q(2))*(dq(1) + dq(2));
    coriolis(2,0) = -dq(0)*(I3_xy*cos(2*q(1) + 2*q(2)) + (I3_xx*sin(2*q(1) + 2*q(2)))/2 - (I3_yy*sin(2*q(1) + 2*q(2)))/2 + (pow(L3_ycom,2)*m3*sin(2*q(1) + 2*q(2)))/2 - (L23_y*L3_ycom*m3*sin(q(2)))/2 + (L23_y*L3_ycom*m3*sin(2*q(1) + q(2)))/2);
    coriolis(2,1) = L23_y*L3_ycom*m3*dq(1)*sin(q(2));
    coriolis(2,2) = 0;

    return coriolis;
}

Eigen::Vector3d DynamicsSolver::Gravity(const Eigen::Vector3d &q)
{
    Eigen::Vector3d gravity;

    gravity(0) = 0;
    gravity(1) = -g*m3*(L3_ycom*sin(q(1) + q(2)) + L23_y*sin(q(1))) - L2_ycom*g*m2*sin(q(1));
    gravity(2) = -L3_ycom*g*m3*sin(q(1) + q(2));

    return gravity;

}