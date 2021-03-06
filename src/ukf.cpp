#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

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
  P_.setIdentity();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .6;
  
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
    
  is_initialized_ = false;
  // dimension of state vector
  n_x_ = 5;
  // dimension of augmented state vector
  n_aug_ = 7;
  // parameter
  lambda_ = 3 - n_aug_;
  // sigma points prediction matrix
  Xsig_pred_ = MatrixXd(n_x_, n_aug_ * 2 + 1);
  // time
  time_us_ = 0.;

  // set weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1 ; i++) weights_(i) = 0.5 / (n_aug_ + lambda_);
  
  //ladar measurement noise covariance matrix
  R_ladar_ = MatrixXd(2, 2);
  R_ladar_ << std_laspx_* std_laspx_, 0,
              0, std_laspy_ * std_laspy_ ;
  
  // radar measurement noise covariance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_* std_radr_, 0, 0,
             0, std_radphi_ * std_radphi_, 0,
             0, 0, std_radrd_*std_radrd_;
  
  // Augmented state covariance matrix
  Q_ = MatrixXd(n_aug_, n_aug_);
  Q_.fill(0.0);
  Q_(5,5) = std_a_ * std_a_;
  Q_(6,6) = std_yawdd_ * std_yawdd_;
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
  if (!is_initialized_) {
    // first measurement:
    cout << "initialize" << endl;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho, phi;
      rho = meas_package.raw_measurements_[0];
      phi = meas_package.raw_measurements_[1];
      x_ << rho * cos(phi), rho * sin(phi), 0, 0, 0;
      cout <<"radar" << endl;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            0, 0, 0;
      cout <<"laser" << x_ << endl;
    }

    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }
  
  //cout << "debug 0" << endl;
  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;
  //cout << "debug 1 " << delta_t <<  endl;

  if (use_radar_==true && 
      meas_package.sensor_type_ == MeasurementPackage::RADAR) { //convert the prediction to radar space
    Prediction(delta_t);
    UpdateRadar(meas_package);
    return;
  }
  
  else if (use_laser_==true &&
           meas_package.sensor_type_ == MeasurementPackage::LASER) {
    Prediction(delta_t);
    UpdateLidar(meas_package);
    return;
  }

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
  
  // augmented sigma points generation
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  VectorXd x_aug = VectorXd(n_aug_);
 
  // first 5 elements of x_aug is (px, py, v, yaw, yawdot)
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  Q_.topLeftCorner(n_x_, n_x_) = P_;
    
  MatrixXd L = Q_.llt().matrixL();
    
  Xsig_aug.col(0) = x_aug;
    for (int i = 0; i < n_aug_; i++) {
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug -sqrt(lambda_ + n_aug_) * L.col(i);
  }
  
  // cout << "Xsig_aug " << endl << Xsig_aug << endl;

  // predict sigma points  
  VectorXd a = VectorXd(n_x_); // first additive term, from integral
  VectorXd b = VectorXd(n_x_); // second additive term, from noise
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {

    const double p_x = Xsig_aug(0,i);
    const double p_y = Xsig_aug(1,i);
    const double v = Xsig_aug(2,i);
    const double yaw = Xsig_aug(3,i);
    const double yawd = Xsig_aug(4,i);
    const double nu_a = Xsig_aug(5,i);
    const double nu_yawdd = Xsig_aug(6,i);
    
    double px_p, py_p;
    
    
    if (fabs(yawd) > 0.001) {  // if yaw_dot != 0
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }
    
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;
    
    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;
    
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
    
    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p; 
  }
  //cout << " Xsig_pred_ " << endl << Xsig_pred_ << endl;
  // predict mean and covariance using predicted sigma points
  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //iterate over sigma points
    x_ += weights_(i) * Xsig_pred_.col(i);
  }
    
  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {      
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
    P_ += weights_(i) * x_diff * x_diff.transpose();
    }
  //cout << "prediction" << endl <<  x_ << endl;
  //cout << "P_ " << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Given prediction and Lidar measurement, update state vector x_ and covariance P_.
  calculate the lidar NIS.
  */
  
  VectorXd z_pred = VectorXd(2);
  z_pred.fill(0.0);
  MatrixXd S = MatrixXd(2,2);
  S.fill(0.0);
  MatrixXd Zsig = MatrixXd(2, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  Zsig = Xsig_pred_.topLeftCorner(2, 2 * n_aug_ + 1);
  // cout << "Xsig_pred_ " << endl  << Xsig_pred_ << endl;
  // cout << "Zsig" << endl << Zsig << endl;
  
  //mean prediction
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  // covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    S += weights_(i) * z_diff * z_diff.transpose();
    
  }

  S = S + R_ladar_;
  
  // Update Step
  VectorXd z = VectorXd(2);
  z << meas_package.raw_measurements_[0], 
       meas_package.raw_measurements_[1];

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, 2);

  //calculate cross correlation btw sigma points in
  // state space and measurement space
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

  //residual
  VectorXd z_diff = Zsig.col(i) - z_pred;

  // state difference
  VectorXd x_diff = Xsig_pred_.col(i) - x_;

  Tc +=  weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, 2);
  MatrixXd S_inv = S.inverse();
  K = Tc * S_inv;

  // update mean and Covariance Matrix
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();
  
  double NIS = (z - z_pred).transpose() * S_inv * (z - z_pred); 
  cout << "NIS " << NIS << endl;

}

MatrixXd UKF::SigmaPointsToRadarSpace(int n_z) {
  /**
   Convert sigma points predicted by Prediction() function to
   Radar Space. Returns a Matrix Zsig
   **/
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double  v = Xsig_pred_(2,i);
    double psi= Xsig_pred_(3,i);
    double v1 = cos(psi)*v;
    double v2 = sin(psi)*v;
    
    double rho = sqrt(px*px + py*py);// if r sufficiently large:
    if (rho > 0.01) {
      Zsig(0,i) = rho;
      Zsig(1,i) = atan2(py, px);
      Zsig(2,i) = (px*v1 + py*v2)/rho;
    }
    else {
      Zsig(0,i) = rho;
      Zsig(1,i) = 0.;
      Zsig(2,i) = v;  // has to be fixed
    }
  }
  return Zsig;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

   Given prediction and Radar measurement, update state vector x_ and covariance P_.
   calculate the radar NIS.
  */
  
  // Prediction in Radar Space
  VectorXd z_pred = VectorXd(3);
  z_pred.fill(0.0);
  MatrixXd S = MatrixXd(3,3);
  S.fill(0.0);
  //create matrix for sigma points in measurement space
  int n_z = 3;
  MatrixXd Zsig = SigmaPointsToRadarSpace(n_z);

  // mean predicted measuremet using Zsig
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
      //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;
      
      //angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
      
      S += weights_(i) * z_diff * z_diff.transpose();
      //cout << "S" << "i" << i << endl << S << endl;
  }
  
  //add measurement noise covariance matrix

  S = S + R_radar_;
  //cout << "S" << endl << S << endl;
  //cout << "R" << endl << R << endl;
  //cout << "z_pred" << endl << z_pred << endl;
  
  // Update Step
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], 
       meas_package.raw_measurements_[1],
       meas_package.raw_measurements_[2];

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, 3);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc +=  weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  MatrixXd S_inv = S.inverse();
  K = Tc * S_inv;

  // update mean and Covariance Matrix
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();
  
  double NIS = (z - z_pred).transpose() * S_inv * (z - z_pred); 
  cout << "NIS " << NIS << endl;

}
