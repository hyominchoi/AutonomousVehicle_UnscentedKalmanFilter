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
  Eigen::VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  if (estimations.size() != ground_truth.size() ||
      estimations.size() ==0) {
    cout << "Invalid estimation or ground truth data" << endl;
    return rmse;
    }
    
  for (int i = 0; i < estimations.size(); ++ i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    VectorXd residual_2 = residual.array() * residual.array();
    rmse += residual_2;
    }
  //cout << rmse << "debug 2" << endl;
  rmse = rmse / estimations.size();
  //cout << rmse << "debug 3" << endl;
  rmse = rmse.array().sqrt();
  //cout <<rmse << "debug 4" << endl;
  return rmse;
}
