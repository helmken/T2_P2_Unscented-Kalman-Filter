#include <iostream>
#include <numeric>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

using namespace std;


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
    VectorXd rmse = VectorXd::Zero(4);

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
    rmse = rmse / estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

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


double NormalizeAngle(double  phi)
{
    return atan2(sin(phi), cos(phi));
}

VectorXd ReadGroundTruth(istringstream& iss)
{
    const int groundTruthDim(6);
    VectorXd groundTruth(groundTruthDim);
    for (int i(0); i < groundTruthDim; ++i)
    {
        iss >> groundTruth[i];
    }

    return groundTruth;
}

VectorXd CalcStandardDeviation(const std::vector<VectorXd> measurements)
{
    VectorXd zeroVec;

    if (measurements.empty())
    {
        printf("CalcStandardDeviation: measurements vector is empty\n");
        return zeroVec;
    }

    const VectorXd firstVec = measurements[0];
    const int measurementDim = firstVec.rows();

    zeroVec = VectorXd::Zero(measurementDim);
    
    const double measurementsCount(measurements.size());

    const VectorXd measurementsSum = std::accumulate(measurements.begin(), measurements.end(), zeroVec);
    const VectorXd measurementsMean = measurementsSum / measurementsCount;

    std::vector<VectorXd> measurementsSquaredDiff;
    for (VectorXd x : measurements)
    {
        VectorXd diff = x - measurementsMean;
        VectorXd diffSquared = diff.array().pow(2);
        measurementsSquaredDiff.push_back(diffSquared);
    }

    const VectorXd measurementsSquaredDiffSum = std::accumulate(
        measurementsSquaredDiff.begin(), 
        measurementsSquaredDiff.end(),
        zeroVec);

    const VectorXd variance = measurementsSquaredDiffSum / measurementsCount;
    const VectorXd stdDev = variance.array().sqrt();

    return stdDev;
}

void CalcStandardDeviation(
    const std::vector<VectorXd> groundTruthValues, 
    double& stdAcc, double& minAcc, double& maxAcc,
    double& stdYawRate, double& minYawRate, double& maxYawRate)
{
    // ground truth : x, y, vx, vy, yaw, yawrate

    minAcc = std::numeric_limits<double>::max();
    maxAcc = std::numeric_limits<double>::min();

    minYawRate = std::numeric_limits<double>::max();
    maxYawRate = std::numeric_limits<double>::min();

    if (groundTruthValues.empty())
    {
        printf("CalcStandardDeviation given ground truth values are empty\n");
        return;
    }

    const double microSecs2Secs = 0.000001;

    // longAcc.size() will be yawRates.size() - 1 because longitudinal accelerations
    // have to be calculated based on deltaT whereas yawRates directly measured
    vector<double> yawRates;
    vector<double> longAccs;

    // ground truth vector = [timestamp, x, y, vx, vy, yaw, yawrate]
    const VectorXd& firstGroundTruth = groundTruthValues[0];
    yawRates.push_back(firstGroundTruth[6]);
    double lastTimeStamp = firstGroundTruth[0] * microSecs2Secs;
    double lastVelocity = sqrt(pow(firstGroundTruth[3], 2) + pow(firstGroundTruth[4], 2));

    for (unsigned int groundTruthIdx(1); groundTruthIdx < groundTruthValues.size(); ++groundTruthIdx)
    {
        const VectorXd& groundTruth = groundTruthValues[groundTruthIdx];
        double curTimeStamp = groundTruth[0] * microSecs2Secs;
        double deltaT = curTimeStamp - lastTimeStamp;

        if (deltaT < 0.001)
        {
            printf("skipping measurement, deltaT=%.f < 0.001\n", deltaT);
        }
        else
        {
            double curVelocity = sqrt(pow(groundTruth[3], 2) + pow(groundTruth[4], 2));
            double curAcceleration = (curVelocity - lastVelocity) / deltaT;
            
            if (curAcceleration < minAcc)
            {
                minAcc = curAcceleration;
            }
            if (curAcceleration > maxAcc)
            {
                maxAcc = curAcceleration;
            }

            printf("curVelocity=%f, lastVelocity= %f, deltaT=%f, curAcceleration=%f\n", 
                curVelocity, lastVelocity, deltaT, curAcceleration);

            longAccs.push_back(curAcceleration);

            yawRates.push_back(groundTruth[6]);

            if (groundTruth[6] < minYawRate)
            {
                minYawRate = groundTruth[6];
            }
            if (groundTruth[6] > maxYawRate)
            {
                maxYawRate = groundTruth[6];
            }

            lastVelocity = curVelocity;
            lastTimeStamp = curTimeStamp;
        }
    }

    const double longAccsSum = std::accumulate(longAccs.begin(), longAccs.end(), 0.0);
    const double longAccsMean = longAccsSum / longAccs.size();

    std::vector<double> longAccsSquaredDiff(longAccs.size());
    std::transform(
        longAccs.begin(), longAccs.end(),
        longAccsSquaredDiff.begin(),
        [longAccsMean](double x)
        {
            double squaredDiff = pow(x - longAccsMean, 2);
            return squaredDiff;
        });

    const double longAccsSquaredDiffSum = std::accumulate(
        longAccsSquaredDiff.begin(), longAccsSquaredDiff.end(), 0.0);

    stdAcc = sqrt(longAccsSquaredDiffSum / longAccsSquaredDiff.size());

    const double yawRatesSum = std::accumulate(yawRates.begin(), yawRates.end(), 0.0);
    const double yawRateMean = yawRatesSum / yawRates.size();

    std::vector<double> yawRatesSquaredDiff(yawRates.size());
    std::transform(
        yawRates.begin(), yawRates.end(),
        yawRatesSquaredDiff.begin(),
        [yawRateMean](double x)
        {
            double squaredDiff = pow(x - yawRateMean, 2);
            return squaredDiff;
        });

    const double yawRatesSquaredDiffSum = std::accumulate(
        yawRatesSquaredDiff.begin(), yawRatesSquaredDiff.end(), 0.0);

    stdYawRate = sqrt(yawRatesSquaredDiffSum / yawRatesSquaredDiff.size());
}
