#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


UKF::UKF()
    :
    is_initialized_(false),
    use_laser_(true),
    use_radar_(true),
    previous_timestamp_(0)
{
    // initial state vector
    x_ = VectorXd::Zero(x_dim);

    // initial covariance matrix
    P_ = MatrixXd::Zero(x_dim, x_dim);

    //create vector for weights
    weights_ = VectorXd(num_sigma_points);

    // set weights
    weights_.fill(0.5 / (lambda + x_aug_dim));
    weights_(0) = lambda / (lambda + x_aug_dim);

    // measurement noise covariance for radar measurements
    R_radar_ = MatrixXd::Zero(radarMeas_dim, radarMeas_dim);
    R_radar_(0, 0) = pow(std_radar_rho, 2);
    R_radar_(1, 1) = pow(std_radar_phi, 2);
    R_radar_(2, 2) = pow(std_radar_rhodot, 2);

    // measurement noise covariance for laser measurements
    R_laser_ = MatrixXd::Zero(laserMeas_dim, laserMeas_dim);
    R_laser_(0, 0) = pow(std_laser_px, 2);
    R_laser_(1, 1) = pow(std_laser_py, 2);

    // map x and y of state into laser measurement space
    H_laser_ = MatrixXd::Zero(laserMeas_dim, x_dim);
    H_laser_(0, 0) = 1.0;
    H_laser_(1, 1) = 1.0;
}

UKF::~UKF()
{
}

void UKF::Initialize(const MeasurementPackage &measurement_pack)
{
    // first measurement
    x_.fill(0.0);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
        // Convert radar from polar to cartesian coordinates and initialize state.
        const VectorXd cartesianPos = Polar2Cartesian(
            measurement_pack.raw_measurements_[0],
            measurement_pack.raw_measurements_[1],
            measurement_pack.raw_measurements_[2]);

        printf("initial measurement R: (%.4f,%.4f)\n", cartesianPos[0], cartesianPos[1]);

        // x_: [pos1 pos2 vel_abs yaw_angle yaw_rate]
        x_[0] = cartesianPos[0]; // p_x
        x_[1] = cartesianPos[1]; // p_y
        //   sqrt(cos ^ 2(phi)* ro_dot ^ 2 + sin ^ 2(phi)* ro_dot ^ 2) 
        // = sqrt(ro_dot ^ 2 * (cos ^ 2(phi)+sin ^ 2(phi))
        // = sqrt(ro_dot ^ 2) 
        // = abs(ro_dot)        
        x_[2] = abs(measurement_pack.raw_measurements_[2]);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
        // Initialize state.
        printf("initial measurement L: (%.4f,%.4f)\n",
            measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1]);

        x_[0] = measurement_pack.raw_measurements_[0];
        x_[1] = measurement_pack.raw_measurements_[1];
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    P_ = MatrixXd::Identity(x_dim, x_dim);

    // done initializing, no need to predict or update
    is_initialized_ = true;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
    if (!is_initialized_)
    {
        Initialize(meas_package);
        return;
    }

    /*****************************************************************************
    * Prediction: predict sigma points, the state, and the state covariance matrix
    ****************************************************************************/

    const double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
    if (fabs(dt) < std::numeric_limits<double>::epsilon())
    {
        printf("***\n*** dt is zero -> no prediction necessary ***\n***\n");
        return;
    }

    MatrixXd augmentedSigmaPoints = GenerateAugmentedSigmaPoints(x_, P_);
    MatrixXd predictedSigmaPoints = CalculatePredictedSigmaPoints(augmentedSigmaPoints, dt);

    VectorXd x_pred; // predicted mean
    MatrixXd P_pred; // predicted covariance
    PredictMeanAndCovariance(predictedSigmaPoints, x_pred, P_pred);

    /*****************************************************************************
    *  Update
    ****************************************************************************/
    // Use the sensor type to perform the update step, update the state and covariance 
    // matrices.
    if (   MeasurementPackage::RADAR == meas_package.sensor_type_
        && use_radar_)
    {
        VectorXd z = VectorXd(radarMeas_dim); // current measurement
        z[0] = meas_package.raw_measurements_[0]; // rho
        z[1] = NormalizeAngle(meas_package.raw_measurements_[1]); // phi
        z[2] = meas_package.raw_measurements_[2]; // rho_dot

        MatrixXd predSigmaPointsInMeasSpace; // predicted sigma points in measurement space
        VectorXd z_pred; // predicted measurement mean
        MatrixXd S_pred; // predicted measurement variance
        PredictRadarMeasurement(
            predictedSigmaPoints, predSigmaPointsInMeasSpace, 
            z_pred, S_pred);

        UpdateStateRadar(
            predictedSigmaPoints, x_pred, P_pred,
            predSigmaPointsInMeasSpace, z_pred, S_pred,
            z, x_, P_);

        radarNisValues_.push_back(CalculateRadarNIS(z, z_pred, S_pred));
        CalculateNisConsistencyFromMeasurements(radar_nis_threshold, radarNisValues_);

        previous_timestamp_ = meas_package.timestamp_;
    }
    else if (   MeasurementPackage::LASER == meas_package.sensor_type_
             && use_laser_)
    {
        VectorXd z = VectorXd(laserMeas_dim);
        z[0] = meas_package.raw_measurements_[0]; // p_x
        z[1] = meas_package.raw_measurements_[1]; // p_y 

        UpdateStateLaser(
            x_pred, P_pred, z, 
            H_laser_, R_laser_,
            x_, P_);

        previous_timestamp_ = meas_package.timestamp_;
    }
}

MatrixXd GenerateAugmentedSigmaPoints(const VectorXd& x, const MatrixXd& P)
{
    // create augmented mean state
    VectorXd x_aug = VectorXd::Zero(x_aug_dim);
    x_aug.head(x_dim) = x;

    // create augmented state covariance
    MatrixXd P_aug = MatrixXd::Zero(x_aug_dim, x_aug_dim);
    P_aug.topLeftCorner(x_dim, x_dim) = P;
    P_aug(x_dim, x_dim) = std_a * std_a;
    P_aug(x_dim + 1, x_dim + 1) = std_yaw_dot * std_yaw_dot;

    // calculate square root matrix
    MatrixXd sqrt_P_aug = P_aug.llt().matrixL();

    // create matrix for augmented sigma points
    MatrixXd sqrt_lambda_Paug = sqrt(lambda + x_aug_dim) * sqrt_P_aug;

    // create augmented sigma point matrix
    MatrixXd augmentedSigmaPoints = MatrixXd(x_aug_dim, num_sigma_points);
    augmentedSigmaPoints.col(0) = x_aug;
    for (int i(0); i < x_aug_dim; ++i)
    {
        augmentedSigmaPoints.col(i + 1) =             x_aug + sqrt_lambda_Paug.col(i);
        augmentedSigmaPoints.col(i + 1 + x_aug_dim) = x_aug - sqrt_lambda_Paug.col(i);
    }

    return augmentedSigmaPoints;
}

MatrixXd CalculatePredictedSigmaPoints(const MatrixXd& augmentedSigmaPoints, const double delta_t)
{
    //create matrix with predicted sigma points as columns
    MatrixXd predictedSigmaPoints = MatrixXd(x_dim, num_sigma_points);

    //predict sigma points
    for (int i(0); i < num_sigma_points; ++i)
    {
        const VectorXd sigmaPt_org = augmentedSigmaPoints.col(i);
        const double px = sigmaPt_org[0];
        const double py = sigmaPt_org[1];
        const double v = sigmaPt_org[2];
        const double psi = sigmaPt_org[3];
        const double psi_dot = sigmaPt_org[4];
        const double nu_a = sigmaPt_org[5];
        const double nu_psi_dot = sigmaPt_org[6];

        VectorXd sigmaPt_pred = VectorXd(x_dim);
        if (fabs(psi_dot) < 0.001)
        {
            sigmaPt_pred[0] = px
                            + v * cos(psi) * delta_t
                            + 0.5 * pow(delta_t, 2) * cos(psi) * nu_a;
            sigmaPt_pred[1] = py
                            + v * sin(psi) * delta_t
                            + 0.5 * pow(delta_t, 2) * sin(psi) * nu_a;
            sigmaPt_pred[2] = v + delta_t * nu_a;
            sigmaPt_pred[3] = psi + 0.5 * pow(delta_t, 2) * nu_psi_dot;
            sigmaPt_pred[4] = psi_dot + delta_t * nu_psi_dot;
        }
        else
        {
            sigmaPt_pred[0] = px
                            + (v / psi_dot) * (sin(psi + psi_dot * delta_t) - sin(psi))
                            + 0.5 * pow(delta_t, 2) * cos(psi) * nu_a;
            sigmaPt_pred[1] = py
                            + (v / psi_dot) * (-cos(psi + psi_dot * delta_t) + cos(psi))
                            + 0.5 * pow(delta_t, 2) * sin(psi) * nu_a;
            sigmaPt_pred[2] = v + delta_t * nu_a;
            sigmaPt_pred[3] = psi
                            + psi_dot * delta_t
                            + 0.5 * pow(delta_t, 2) * nu_psi_dot;
            sigmaPt_pred[4] = psi_dot + delta_t * nu_psi_dot;
        }

        predictedSigmaPoints.col(i) = sigmaPt_pred;
    }

    return predictedSigmaPoints;
}

void UKF::PredictMeanAndCovariance(
    const MatrixXd& predictedSigmaPoints, VectorXd& x_pred, MatrixXd& P_pred)
{
    //predicted state mean
    x_pred = VectorXd::Zero(x_dim);
    for (int i(0); i < num_sigma_points; ++i)
    {  
        //iterate over sigma points
        x_pred += weights_(i) * predictedSigmaPoints.col(i);
    }

    //predicted state covariance matrix
    P_pred = MatrixXd::Zero(x_dim, x_dim);
    for (int i(0); i < num_sigma_points; ++i)
    {  
        // state difference
        VectorXd x_diff = predictedSigmaPoints.col(i) - x_pred;
        
        //angle normalization
        x_diff(3) = NormalizeAngle(x_diff(3));
        
        P_pred += weights_(i) * x_diff * x_diff.transpose();
    }
}

void UKF::PredictRadarMeasurement(
    const MatrixXd& predSigmaPoints,
    MatrixXd& predSigmaPointsInMeasSpace,
    VectorXd& z_pred,
    MatrixXd& S_pred) 
{
    //create matrix for sigma points in measurement space
    predSigmaPointsInMeasSpace = MatrixXd::Zero(radarMeas_dim, num_sigma_points);

    //transform sigma points into measurement space
    for (int i(0); i < num_sigma_points; ++i)
    {
        const VectorXd& predSigmaPt = predSigmaPoints.col(i);
        
        const double px = predSigmaPt[0];
        const double py = predSigmaPt[1];
        const double v = predSigmaPt[2];
        const double psi = predSigmaPt[3];
        const double rho = sqrt(px * px + py * py);
        
        double phi(0.0);
        if (fabs(phi) > 0.00001)
        {
            phi = atan2(py, px);
        }

        double rhoDot(0.0);
        if (fabs(rho) > 0.00001)
        {
            rhoDot =   (px * cos(psi) * v + py * sin(psi) * v)
                     / rho;
        }

        predSigmaPointsInMeasSpace.col(i)[0] = rho;
        predSigmaPointsInMeasSpace.col(i)[1] = phi;
        predSigmaPointsInMeasSpace.col(i)[2] = rhoDot;
    }

    //mean predicted measurement
    z_pred = VectorXd::Zero(radarMeas_dim);

    //calculate mean predicted measurement
    for (int i(0); i < num_sigma_points; ++i)
    {
        z_pred = z_pred + weights_[i] * predSigmaPointsInMeasSpace.col(i);
    }

    // normalize angle phi
    z_pred[1] = NormalizeAngle(z_pred[1]);

    //measurement covariance matrix S
    S_pred = MatrixXd::Zero(radarMeas_dim, radarMeas_dim);

    //calculate measurement covariance matrix S
    for (int i(0); i < num_sigma_points; ++i)
    {
        VectorXd diff = predSigmaPointsInMeasSpace.col(i) - z_pred;
        S_pred = S_pred + weights_[i] * (diff * diff.transpose());
    }

    // add radar measurement noise R
    S_pred = S_pred + R_radar_;
}

void UKF::UpdateStateRadar(
    const MatrixXd& predictedSigmaPoints, VectorXd& x_pred, MatrixXd& P_pred,
    MatrixXd& predSigmaPointsInMeasSpace, VectorXd& z_pred, MatrixXd& S_pred,
    VectorXd& z,
    VectorXd& x_updated, MatrixXd& P_updated)
{
    // Use radar data to update the belief about the object's position. 

    //calculate cross correlation matrix Tc
    MatrixXd Tc = MatrixXd::Zero(x_dim, radarMeas_dim);
    for (int i(0); i < num_sigma_points; ++i)
    {
        // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
        VectorXd xDiff = predictedSigmaPoints.col(i) - x_pred;
        VectorXd zDiff = predSigmaPointsInMeasSpace.col(i) - z_pred;
        Tc = Tc + weights_[i] * xDiff * (zDiff.transpose());
    }

    //calculate Kalman gain K;
    MatrixXd K = Tc * (S_pred.inverse());

    //update state mean and covariance matrix
    VectorXd measDiff = z - z_pred;
    measDiff[1] = NormalizeAngle(measDiff[1]);
    
    x_updated = x_pred + K * measDiff;
    P_updated = P_pred - K * S_pred * (K.transpose());
}

void UKF::PredictLaserMeasurement(
    const MatrixXd& predSigmaPoints,
    MatrixXd& predSigmaPointsInMeasSpace,
    VectorXd& z_pred,
    MatrixXd& S_pred)
{
    //create matrix for sigma points in measurement space
    predSigmaPointsInMeasSpace = MatrixXd::Zero(laserMeas_dim, num_sigma_points);

    //transform sigma points into measurement space
    for (int i(0); i < num_sigma_points; ++i)
    {
        const VectorXd& predSigmaPt = predSigmaPoints.col(i);

        const double px = predSigmaPt[0];
        const double py = predSigmaPt[1];

        predSigmaPointsInMeasSpace.col(i)[0] = px;
        predSigmaPointsInMeasSpace.col(i)[1] = py;
    }

    //mean predicted measurement
    z_pred = VectorXd::Zero(laserMeas_dim);

    //calculate mean predicted measurement
    for (int i(0); i < num_sigma_points; ++i)
    {
        z_pred = z_pred + weights_[i] * predSigmaPointsInMeasSpace.col(i);
    }

    //measurement covariance matrix S
    S_pred = MatrixXd::Zero(laserMeas_dim, laserMeas_dim);

    //calculate measurement covariance matrix S
    for (int i(0); i < num_sigma_points; ++i)
    {
        VectorXd diff = predSigmaPointsInMeasSpace.col(i) - z_pred;
        S_pred = S_pred + weights_[i] * (diff * diff.transpose());
    }

    // add laser measurement noise R
    S_pred = S_pred + R_laser_;
}

void UKF::UpdateStateLaser(
    const VectorXd& x_pred, const MatrixXd& P_pred,
    const VectorXd& z,
    const MatrixXd& H, const MatrixXd& R,
    VectorXd& x_updated, MatrixXd& P_updated)
{
    const VectorXd z_pred = H * x_pred;
    const VectorXd y = z - z_pred;

    const MatrixXd Ht = H.transpose();
    const MatrixXd S = H * P_pred * Ht + R;
    const MatrixXd Si = S.inverse();
    const MatrixXd PHt = P_pred * Ht;
    const MatrixXd K = PHt * Si;

    //new estimate
    x_updated = x_pred + (K * y);
    const long x_size = x_pred.size();
    const MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_updated = (I - K * H) * P_pred;
}


double CalculateRadarNIS(const VectorXd& z, const VectorXd& z_pred, const MatrixXd& S_pred)
{
    VectorXd zDiff = z - z_pred;
    double nis = zDiff.transpose() * S_pred.inverse() * zDiff;
    return nis;
}

double CalculateLaserNIS(const VectorXd& z, const VectorXd& z_pred, const MatrixXd& S_pred)
{
    VectorXd zDiff = z - z_pred;
    double nis = zDiff.transpose() * S_pred.inverse() * zDiff;
    return nis;
}

void UKF::CalculateNisConsistency()
{
    printf("radar NIS:\n");
    CalculateNisConsistencyFromMeasurements(radar_nis_threshold, radarNisValues_);

    printf("laser NIS:\n");
    CalculateNisConsistencyFromMeasurements(laser_nis_threshold, laserNisValues_);
}

double CalculateNisConsistencyFromMeasurements(const double nisThreshold, const std::vector<double>& nisValues)
{
    if (0 == nisValues.size())
    {
        printf("given nisValues are empty -> no calculation possible");
        return 0.0;
    }

    int aboveThreshold(0);
    for (const double& nisValue : nisValues)
    {
        if (nisValue > nisThreshold)
        {
            ++aboveThreshold;
        }
    }
    const double fraction = (double)aboveThreshold / (double)nisValues.size();

    //printf("CalculateNisConsistency: %i of %lu measurements are above %.3f -> fraction=%.3f\n",
    //    aboveThreshold, nisValues.size(), nisThreshold, fraction);

    return fraction;
}

