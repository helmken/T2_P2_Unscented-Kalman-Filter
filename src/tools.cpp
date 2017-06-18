#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools()
{
}

Tools::~Tools()
{
}

VectorXd Tools::CalculateRMSE(
    const vector<VectorXd> &estimations,
    const vector<VectorXd> &ground_truth)
{
    // TODO: Calculate the RMSE here. (DONE)

    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    // ... your code here
    if (estimations.size() != ground_truth.size()
        || estimations.empty())
    {
        cout << "CalculateRMSE: invalid sizes of estimations and ground truths" << endl;
        return rmse;
    }

    //accumulate squared residuals
    for (size_t i = 0; i < estimations.size(); ++i)
    {
        // ... your code here
        const VectorXd &curEstimate = estimations[i];
        const VectorXd &curGroundTruth = ground_truth[i];

        if (curEstimate.size() != curGroundTruth.size()
            || 0 == curEstimate.size())
        {
            cout << "CalculateRMSE: invalid size of estimation or ground truth" << endl;
            return rmse;
        }

        VectorXd residual = curEstimate - curGroundTruth;
        VectorXd residual_sq = residual.array() * residual.array();

        rmse += residual_sq;
    }

    //calculate the mean
    // ... your code here
    rmse = rmse / estimations.size();

    //calculate the squared root
    // ... your code here
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

Eigen::VectorXd Polar2Cartesian(
    const double rho, const double phi, const double rho_dot)
{
    VectorXd result = VectorXd(4);

    result[0] = rho * cos(phi); // p_x
    result[1] = rho * sin(phi); // p_y
    result[2] = rho_dot * cos(phi); // approximation of v_x
    result[3] = rho_dot * sin(phi); // approximation of v_y

    return result;
}

double NormalizeAngle(double radian)
{
    if (radian >= -M_PI && radian <= M_PI)
    {
        return radian;
    }

    double normalized(radian);
    while (normalized > M_PI)
    {
        normalized -= 2 * M_PI;
    }

    while (normalized < -M_PI)
    {
        normalized += 2 * M_PI;
    }

    //printf("normalizing angle %.5f to %.5f\n", radian * rad2deg, normalized * rad2deg);
    
    return normalized;
}
