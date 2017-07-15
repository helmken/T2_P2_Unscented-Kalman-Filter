#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>


using Eigen::MatrixXd;
using Eigen::VectorXd;


// State dimension
const int x_dim(5);

// Augmented state dimension
const int x_aug_dim(7);

// Sigma point spreading parameter
const double lambda(3 - x_aug_dim);

const int num_sigma_points(2 * x_aug_dim + 1);

// set radar measurement dimension: radar can measure r, phi, and r_dot
const int radarMeas_dim(3);

// set laser measurement dimension: laser measures position x and y
const int laserMeas_dim(2);

// Process noise standard deviation longitudinal acceleration in m/s^2
const double std_a(0.5);

// Process noise standard deviation yaw acceleration in rad/s^2
const double std_yaw_dot(M_PI / 4.0);

// rubric values: [.09, .10, .40, .30].

// Laser measurement noise standard deviation position1 in m
const double std_laser_px(0.15);

// Laser measurement noise standard deviation position2 in m
const double std_laser_py(0.15);

// Radar measurement noise standard deviation radius in m
const double std_radar_rho(0.3);    // value from lesson: 0.3

// Radar measurement noise standard deviation angle in rad
const double std_radar_phi(0.03);   // value from lesson: 0.0175

// Radar measurement noise standard deviation radius change in m/s
const double std_radar_rhodot(0.3); // value from lesson: 0.1

// radar measurements are 3D, so NIS consistency threshold is 7.815
const double radar_nis_threshold(7.815);

// laser measurements are 2D, so NIS consistency threshold is 5.991
const double laser_nis_threshold(5.991);


class UKF
{

private:
    
    int measurement_count_;

    // initially set to false, set to true in first call of ProcessMeasurement
    bool is_initialized_;

    // if this is false, laser measurements will be ignored (except for init)
    bool use_laser_;

    // if this is false, radar measurements will be ignored (except for init)
    bool use_radar_;

    // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    VectorXd x_;

    // state covariance matrix
    MatrixXd P_;

    // measurement covariance noise R for radar measurements
    MatrixXd R_radar_;

    // measurement covariance noise R for laser measurements
    MatrixXd R_laser_;

    // time when the state is true, in us
    long long previous_timestamp_;

    // Weights of sigma points
    VectorXd weights_;

    // track calculated NIS values for laser and radar
    std::vector<double> laserNisValues_;
    std::vector<double> radarNisValues_;

public:
    /**
     * Constructor
     */
    UKF();

    /**
     * Destructor
     */
    virtual ~UKF();

    /**
    * ProcessMeasurement
    * @param meas_package The latest measurement data of either radar or laser
    */
    void ProcessMeasurement(MeasurementPackage meas_package);

    VectorXd GetX()
    {
        return x_;
    };

    /**
     * calculate the percentage of NIS values for laser and radar 
     */
    void CalculateNisConsistency();

private:

    /**
     * Initialization with first measurement
     */
    void Initialize(const MeasurementPackage &measurement_pack);

    void PredictMeanAndCovariance(
        const MatrixXd& predictedSigmaPoints, VectorXd& x_pred, MatrixXd& P_pred);

    void PredictRadarMeasurement(
        const MatrixXd& predSigmaPoints,
        MatrixXd& predSigmaPointsInMeasSpace,
        VectorXd& z_pred,
        MatrixXd& S_pred);

    void UpdateStateRadar(
        const MatrixXd& predictedSigmaPoints, VectorXd& x_pred, MatrixXd& P_pred,
        MatrixXd& predSigmaPointsInMeasSpace, VectorXd& z_pred, MatrixXd& S_pred,
        VectorXd& z,
        VectorXd& x_updated, MatrixXd& P_updated);

    void PredictLaserMeasurement(
        const MatrixXd& predSigmaPoints,
        MatrixXd& predSigmaPointsInMeasSpace,
        VectorXd& z_pred,
        MatrixXd& S_pred);

    void UpdateStateLaser(
        const MatrixXd& predictedSigmaPoints, VectorXd& x_pred, MatrixXd& P_pred,
        MatrixXd& predSigmaPointsInMeasSpace, VectorXd& z_pred, MatrixXd& S_pred,
        VectorXd& z,
        VectorXd& x_updated, MatrixXd& P_updated);
};

MatrixXd GenerateAugmentedSigmaPoints(const VectorXd& x, const MatrixXd& P);

MatrixXd CalculatePredictedSigmaPoints(const MatrixXd& augmentedSigmaPoints, const double delta_t);

double CalculateRadarNIS(const VectorXd& z, const VectorXd& z_pred, const MatrixXd& S_pred);

double CalculateLaserNIS(const VectorXd& z, const VectorXd& z_pred, const MatrixXd& S_pred);

double CalculateNisConsistencyFromMeasurements(const double nisThreshold, const std::vector<double>& nisValues);

#endif /* UKF_H */
