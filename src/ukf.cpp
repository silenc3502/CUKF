#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define EPSILON	0.0001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  /* It comes from part of UKF lecture 14, 17, and 23. */
  n_x_ = x_.size();
  n_aug_ = n_x_ + 2;
  n_sig_ = 2 * n_aug_ + 1;
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  lambda_ = 3 - n_aug;
  weights_ = VectorXd(n_sig_);

  /* Init Measurement Noise Covariance Matrices */
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << pow(std_radr_, 2), 0, 0,
              0, pow(std_radphi_, 2), 0,
              0, 0, pow(std_radrd_, 2);

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << pow(std_laspx_, 2), 0,
              0, pow(std_laspy_, 2);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if(!is_initialized_)
  {
    /* It's come from FusionEKF.cpp and add CTRV Model */
    P_ << 1, 0, 0, 0, 0,
	  0, 1, 0, 0, 0,
	  0, 0, 1, 0, 0,
	  0, 0, 0, 1, 0,
	  0, 0, 0, 0, 1;

    if(measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];

      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);
      float v = sqrt(pow(vx, 2) + pow(vy, 2));

      x_ << px, py, v, 0, 0;
    }
    else if(measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      x_ << measurement_pack.raw_measurements_[0],
            measurement_pack.raw_measurements_[1],
            0,
            0,
            0;

      if(fabs(x_(0)) < EPSILON && fabs(x_(1)) < EPSILON)
      {
        x_(0) = EPSILON;
        x_(1) = EPSILON;
      }
    }

    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for(int i = 1; i < weights_.size(); i++)
      weights_(i) = 0.5 / (n_aug_ + lambda_);

    time_us_ = measurement_pack.timestamp_;
    is_initialized_ = true;

    return;
  }

  double dt = (measurement_pack.timestamp_ - time_us_);
  dt /= 1000000.0;
  time_us_ = measurement_pack.timestamp_;
  Prediction(dt);

  if(measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    UpdateRadar(measurement_pack);

  if(measurement_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    UpdateLidar(measurement_pack);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /* It comes from UKF lecture 17. */
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug, n_sig_);

  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = pow(std_a_, 2);
  P_aug(n_x_ + 1, n_x_ + 1) = pow(std_yawdd_, 2);

  MatrixXd L = P_aug.llt().matrixL();

  XSig_aug.col(0) = x_aug;

  for(int i = 0; i < n_aug; i++)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda + n_aug) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug) = x_aug - sqrt(lambda + n_aug) * L.col(i);
  }

  /* It comes from UKF lecture 20.
     Sigma Point Prediction */
  for(int i = 0; i < n_sig_; i++)
  {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    double px_p, py_p;

    if(fabs(yawd) > EPSILON)
    {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else
    {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    px_p = px_p + 0.5 * nu_a * pow(delta_t, 2) * cos(yaw);
    py_p = py_p + 0.5 * nu_a * pow(delta_t, 2) * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * pow(delta_t, 2);
    yawd_p = yawd_p + nu_yawdd * delta_t;

    Xsig_pred(0, i) = px_p;
    Xsig_pred(1, i) = py_p;
    Xsig_pred(2, i) = v_p;
    Xsig_pred(3, i) = yaw_p;
    Xsig_pred(4, i) = yawd_p;
  }

  /* It comes from UKF lecture 23. */
  x_ = Xsig_pred_ * weigths_;
  P_.fill(0.0);

  for(int i = 0; i < n_sig_; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    while(x_diff(3) > M_PI)
      x_diff(3) -= 2.0 * M_PI;

    while(x_diff(3) < -M_PI)
      x_diff(3) += 2.0 * M_PI;

    P_ = P_ + weights(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  /* It comes from UKF lecture 26 */
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  for(int i = 0; i < n_sig_; i++)
  {
    double p_x = Xsig_pred(0, i);
    double p_y = Xsig_pred(1, i);
    double v = Xsig_pred(2, i);
    double yaw = Xsig_pred(3, i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    Zsig(0, i) = sqrt(pow(p_x, 2) + pow(p_y, 2));
    Zsig(1, i) = atan2(p_y, p_x);
    Zsig(2, i) = (p_x * v1 + p_y * v2) / Zsig(0, i);
  }

  UpdateConditionalUKF(meas_package, Zsig, n_z);
}

/* It comes from UKF lecture 26, and 29. */
void UKF::UpdateConditionalUKF(MeasurementPackage mp, MatrixXd zs, int n_z)
{
  /* Here is 26 part. */
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z, n_z);
  z_pred = zs * weights_;

  S.fill(0.0);

  for(int i = 0; i < n_sig_; i++)
  {
    VectorXd z_diff = zs.col(i) - z_pred;

    while(z_diff(1) > M_PI)
      z_diff(1) -= 2.0 * M_PI;

    while(z_diff(1) < -M_PI)
      z_diff(1) += 2.0 * M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);

  if(mp.sensor_type_ == MeasurementPackage::RADAR)
    R = R_radar_;
  else if(mp.sensor_type_ == MeasurementPackage::LASER)
    R = R_lidar_;

  S = S + R;

  /* Here is 29 part. */
  MatrixXd Tc = MatrixXd(n_x, n_z);
  Tc.fill(0, 0);
  for(int i = 0; i < n_sig_; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if(mp.sensor_type_ == MEASUREMENT::RADAR)
    {
      while(z_diff(1) > M_PI)
        z_diff(1) -= 2.0 * M_PI;

      while(z_diff(1) < -M_PI)
        z_diff(1) += 2.0 * M_PI;
    }

    VectorXd x_diff = Xsig_pred.col(i) - x;

    while(x_diff(3) > M_PI)
      x_diff(3) -= 2.0 * M_PI;

    while(x_diff(3) < -M_PI)
      x_diff(3) += 2.0 * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - z_pred;

  if(mp.sensor_type_ == MeasurementPackage::RADAR)
  {
    while(z_diff(1) > M_PI)
      z_diff(1) -= 2.0 * M_PI;

    while(z_diff(1) < -M_PI)
      z_diff(1) += 2.0 * M_PI;
  }

  x_ = x_ + K * z_diff;
  P_ = P_ - K * s * K.transpose();
}
