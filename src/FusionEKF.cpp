#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  //ekf_.x_ = VectorXd(4); //这里不知道是不是必要的，我感觉不是
  ekf_.P_ = MatrixXd(4, 4); //状态不确定性（相对于测量不确定性R）
  ekf_.P_ << 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 100, 0,
	  0, 0, 0, 100;

  ekf_.F_ = MatrixXd(4,4); //状态转移矩阵
  ekf_.F_ << 1, 0, 0, 0,
			 0, 1, 0, 0,
			 0, 0, 1, 0,
			 0, 0, 0, 1;
  //process noises  Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << 0, 0, 0, 0,
			 0, 0, 0, 0,
			 0, 0, 0, 0,
			 0, 0, 0, 0;
  H_laser_ << 1, 0, 0, 0,
			  0, 1, 0, 0;
 
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

	previous_timestamp_ = measurement_pack.timestamp_;
	
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
		
		float ro;
		float theta;
		float ro_dot;
		
		ro = measurement_pack.raw_measurements_[0];
		theta = measurement_pack.raw_measurements_[1];
		ro_dot = measurement_pack.raw_measurements_[2];
		float px, py, vx, vy;
		px = ro * cos(theta);
		py = ro * sin(theta);
		vx = ro_dot*cos(theta);
		vy = ro_dot*sin(theta);
		
		ekf_.x_ << px, py, vx, vy;

		Hj_ = tools.CalculateJacobian(ekf_.x_);
		
		
		ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);
		
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
		ekf_.x_[0] = measurement_pack.raw_measurements_[0];
		ekf_.x_[1] = measurement_pack.raw_measurements_[1];
		ekf_.x_[2] = 0;
		ekf_.x_[3] = 0;

		H_laser_ << 1, 0, 0, 0,
					0, 1, 0, 0;
		
		ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, H_laser_, R_laser_, ekf_.Q_);
		
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
 
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds

  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt*dt;
  float dt_3 = dt_2*dt;
  float dt_4 = dt_3*dt;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  float noise_ax = 9;
  float noise_ay = 9;
  ekf_.Q_ << dt_4 * noise_ax / 4.0, 0, dt_3*noise_ax / 2.0, 0,
	  0, dt_4*noise_ay / 4.0, 0, dt_3*noise_ay / 2.0,
	  dt_3*noise_ax / 2.0, 0, dt_2*noise_ax, 0,
	  0, dt_3*noise_ay / 2.0, 0, dt_2*noise_ay;

  ekf_.Predict();
  

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  Hj_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.H_ = Hj_;
	  ekf_.R_ = R_radar_;
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);

	  
  } else {
	  
    // Laser updates
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
