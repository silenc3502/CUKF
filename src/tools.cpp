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

  /* It comes from lecture 22 of EKF. */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if(estimations.size() == 0)
  {
    cout << "input is empty" << endl;
    return rmse;
  }

  if(estimations.size() != ground_truth.size())
  {
    cout << "Invalid estimation or ground truth. Data should have the same size" << endl;
    return rmse;
  }

  for(unsigned int i = 0; i < estimations.size(); i++)
  {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}
