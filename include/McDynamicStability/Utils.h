#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <memory>

namespace mc_impact
{

struct ZMPArea
{
  double min_x;
  double max_x;
  double min_y;
  double max_y;
};
struct zmpSupportContact
{
  std::string bodyName;
  std::string sensorName;
};

template<typename Point>
void pointsToInequalityMatrix(const std::vector<Point> & inputPoints,
                              Eigen::MatrixXd & G,
                              Eigen::VectorXd & h,
                              double miniSlope = 0.01,
                              double maxSlope = 1000.0);


template<typename Point>
void pointsToInequalityMatrix(const std::vector<std::shared_ptr<Point>> & inputPoints,
                              Eigen::MatrixXd & G,
                              Eigen::VectorXd & h,
			      std::vector<Eigen::Vector2d> & point,
                              double miniSlope = 0.01,
                              double maxSlope = 1000.0);

Eigen::Matrix3d crossMatrix(const Eigen::Vector3d & input);

} // namespace mc_impact
