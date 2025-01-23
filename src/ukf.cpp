#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 0.25;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.4;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  P_ <<  0.2, 0, 0, 0, 0,
         0, 0.2, 0, 0, 0,
         0, 0, 0.3, 0, 0,
         0, 0, 0, 0.3, 0,
         0, 0, 0, 0, 0.3;

  Q_ = MatrixXd(2,2);
  Q_ << std_a_*std_a_, 0,
       0,     std_yawdd_*std_yawdd_;

  R_radar_ = MatrixXd(3,3);
  R_radar_ <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  R_lidar_ = MatrixXd(2,2);
  R_lidar_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

  is_initialized_ = false;      
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  n_sig_ = 2 * n_aug_ + 1;

  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  weights_ = VectorXd(2*n_aug_ + 1);

  double wi = 0.5/(lambda_ + n_aug_);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for (int i=1; i< n_sig_; i++)
    weights_(i) =  wi;

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_) {
    // first measurement
    std::cout << "UKF: " << std::endl;
    x_ <<  0, 0, 0.5, 0.15, 0.01;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // radar  is converted from polar to cartesian coordinates 
      //         and initialize state
      x_(0) = meas_package.raw_measurements_(0) * cos(meas_package.raw_measurements_(1));
      x_(1) = meas_package.raw_measurements_(0) * sin(meas_package.raw_measurements_(1));
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;  
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  } 
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
    UpdateLidar(meas_package);
  }

  // print the output
  std::cout << "x_ = " << x_ << std::endl;
  std::cout << "P_ = " << P_ << std::endl;

}

MatrixXd UKF::GenerateAugSigmaPoints(){

  /*******************************************************************************
  * Augmented sigma points generation
  ******************************************************************************/

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug.bottomRightCorner(2, 2) = Q_;
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(5) = x_;

  double lplxsqrt = sqrt(lambda_ + n_aug_);
  MatrixXd A = P_aug.llt().matrixL();
  MatrixXd A_scaled = A * lplxsqrt;

  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  Xsig_aug.fill(0.0);
 
  for(int i=0; i< n_sig_; ++i)
  {
    if (i == 0)
      Xsig_aug.col(i) = x_aug;
    else if(i>=1 && i<= n_aug_)
      Xsig_aug.col(i) = x_aug + A_scaled.col((i-1) % n_aug_);
    else
    {
      Xsig_aug.col(i) = x_aug - A_scaled.col((i-1) % n_aug_);  
    }
  }

  return Xsig_aug;
}

MatrixXd UKF::PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t)
{
  /*******************************************************************************
  * Sigma points prediction
  ******************************************************************************/

  MatrixXd Xsig_pred = MatrixXd(n_x_, n_sig_);
  Xsig_pred = Xsig_aug.topRows(n_x_);
  VectorXd noise = VectorXd(n_x_);
  VectorXd integral = VectorXd(n_x_);
  double vel, yaw,yaw_rate, long_acc_noise, yaw_acc_noise, vtoyrt_ratio, delta_yaw;
  for(int i=0; i<n_sig_; i++)
  {
    vel = Xsig_aug(2, i);
    yaw = Xsig_aug(3, i);
    yaw_rate = Xsig_aug(4, i);
    long_acc_noise = Xsig_aug(5, i);
    yaw_acc_noise = Xsig_aug(6, i);

    noise << 
      0.5 * delta_t * delta_t * cos(yaw) * long_acc_noise,
      0.5 * delta_t * delta_t * sin(yaw) * long_acc_noise,
      delta_t * long_acc_noise,
      0.5 * delta_t * delta_t * yaw_acc_noise,
      delta_t * yaw_acc_noise;

    if (fabs(yaw_rate) < 0.001)
    {
      integral << 
        0.5 * vel * cos(yaw) * delta_t,
        0.5 * vel * sin(yaw) * delta_t,
        0,
        yaw_rate * delta_t,
        0;
    }
    else
    {
      vtoyrt_ratio = vel/yaw_rate;
      delta_yaw = yaw + yaw_rate * delta_t;
      integral <<
        vtoyrt_ratio * (sin(delta_yaw) - sin(yaw)),
        vtoyrt_ratio * (cos(yaw) - cos(delta_yaw)),
        0,
        yaw_rate * delta_t,
        0;
    }
    Xsig_pred.col(i) += integral + noise; 
  }
  return Xsig_pred;
}

void UKF::PredictMeanAndCovariance()
{
  /*******************************************************************************
  * Mean and Covariance from predicted Sigma Points
  ******************************************************************************/
  // VectorXd x = VectorXd(n_x_);
  // MatrixXd P = MatrixXd(n_x_, n_x_);
  x_.fill(0.0);
  P_.fill(0.0);
  for(int i=0; i< n_sig_; i++)
  {
    x_ +=  weights_(i) * Xsig_pred_.col(i);
  }

  // predict state covariance matrix
  for(int i=0; i< n_sig_ ; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

}

void UKF::Prediction(double delta_t) {

  MatrixXd Xsig_aug = GenerateAugSigmaPoints();
  Xsig_pred_ = PredictSigmaPoints(Xsig_aug, delta_t);
  PredictMeanAndCovariance();
  
}

std::tuple<VectorXd, MatrixXd , MatrixXd> UKF::PredictLidarMeasurement()
{
  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z,n_z);

  for (int i = 0; i< n_sig_; i++)
  {
    Zsig(0,i) = Xsig_pred_(0,i);  
    Zsig(1,i) = Xsig_pred_(1,i);
  }

  z_pred.fill(0.0);
  for(int i=0; i< n_sig_; i++)
  {
    z_pred +=  weights_(i) * Zsig.col(i);
  }
  
  S.fill(0.0);
  for(int i=0; i< n_sig_; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  S = S + R_lidar_;
  return std::make_tuple(z_pred, S, Zsig);
}


std::tuple<VectorXd, MatrixXd , MatrixXd> UKF::PredictRadarMeasurement()
{
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  for (int i = 0; i< n_sig_; i++)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double psi = Xsig_pred_(3,i);

    double rho = sqrt(px * px + py * py);
    double phi = atan2(py,px);
    double rho_d = (fabs(rho) > 0.001) ? (px * v * cos(psi) + py * v * sin(psi)) / rho : 0; 

    Zsig.col(i) << rho,
                   phi,
                   rho_d;
  }
  
  z_pred.fill(0.0);
  for(int i=0; i< n_sig_; i++)
  {
    z_pred +=  weights_(i) * Zsig.col(i);
  }
  
  S.fill(0.0);
  for(int i=0; i< n_sig_; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R_radar_;
  return std::make_tuple(z_pred, S, Zsig);
}

void UKF::UpdateState(VectorXd z, VectorXd z_pred, MatrixXd S, MatrixXd Zsig, MeasurementPackage::SensorType mtype )
{
  int n_z = (mtype == MeasurementPackage::RADAR) ? 3 : 2;
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  for(int i=0; i< n_sig_; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (mtype == MeasurementPackage::RADAR)
    {
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    }

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - z_pred;
  if (mtype == MeasurementPackage::RADAR) {
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
  }

  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  auto [z_pred,S,Zsig] = PredictLidarMeasurement();
  VectorXd z = VectorXd(2);
  z(0) = meas_package.raw_measurements_(0);
	z(1) = meas_package.raw_measurements_(1);
  UpdateState(z, z_pred, S, Zsig, meas_package.sensor_type_); 
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  auto [z_pred,S, Zsig] = PredictRadarMeasurement();
  VectorXd z = VectorXd(3);
  z(0) = meas_package.raw_measurements_(0);
	z(1) = meas_package.raw_measurements_(1); 	
	z(2) = meas_package.raw_measurements_(2);
  UpdateState(z, z_pred, S, Zsig, meas_package.sensor_type_);
}