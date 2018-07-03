#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd RMSE = VectorXd(4);
	RMSE << 0, 0, 0, 0;
	if (estimations.size() == 0)
	{
		cout << "estimations vector size : 0" << endl;
		return RMSE;
	}
	if (ground_truth.size() != estimations.size())
	{
		cout << "estimations vector size doesn't match ground truth vector size" << endl;
		return RMSE;
	}
	for (int i = 0; i < estimations.size();++i) 
	{
		VectorXd residual = estimations[i] - ground_truth[i];
		VectorXd residual_2 = residual.array()*residual.array();
		RMSE = RMSE + residual_2;
	}
	RMSE = RMSE / estimations.size();
	RMSE = RMSE.array().sqrt();
	return RMSE;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	MatrixXd Hj(3, 4);
	float px, py, vx, vy;
	px = x_state[0];
	py = x_state[1];
	vx = x_state[2];
	vy = x_state[3];

	
	float c1 = px*px + py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);
	
	
	//check division by zero
	if (fabs(c1) < 0.0001) {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		Hj << 0, 0, 0, 0, 
			0, 0, 0, 0,
			0, 0, 0, 0;//这里是否必要？ 
		return Hj;
	}
		
	//compute the Jacobian matrix
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py*(vx*py - vy*px) / c3, px*(px*vy - py*vx) / c3, px / c2, py / c2;

	
	return Hj;
}
