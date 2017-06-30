#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

const double rad2deg(180.0 / M_PI);


class Tools
{
public:
    /**
    * Constructor.
    */
    Tools();

    /**
    * Destructor.
    */
    virtual ~Tools();

    /**
    * A helper method to calculate RMSE.
    */
    VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);
};

/**
* map given 3D polar measurement to 4D cartesian coordinates (px, py, vx, vy)
* @param rho range
* @param phi bearing
* @param rho_dot radial velocity, range rate
* @return 4D vector (px, py, vx, vy)
*/
Eigen::VectorXd Polar2Cartesian(
    const double rho, const double phi, const double rho_dot);

double NormalizeAngle(double radian);

VectorXd ReadGroundTruth(std::istringstream& iss);

VectorXd CalcStandardDeviation(const std::vector<VectorXd> measurements);

void CalcStandardDeviation(const std::vector<VectorXd> groundTruthValues,
    double& stdAcc, double& minAcc, double& maxAcc,
    double& stdYawRate, double& minYawRate, double& maxYawRate);

#endif /* TOOLS_H_ */
