#include "mc_dynamicStability/utils.h"

#include <polytope/staticPoint.h>

template<typename T>
T sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

void clampSlope(double & slope, const double & mini, const double & max)
{
  double sign = sgn<double>(slope);

  if(fabs(slope) <= mini)
  {
    slope = sign * mini;
  }
  else if(fabs(slope) >= max)
  {
    slope = sign * max;
  }
}

template<typename Point>
void mc_impact::pointsToInequalityMatrix(const std::vector<Point> & inputPoints,
                                         Eigen::MatrixXd & G,
                                         Eigen::VectorXd & h,
                                         double miniSlope,
                                         double maxSlope)
{

  int vertexNumber = static_cast<int>(inputPoints.size());
  int dim = static_cast<int>(inputPoints[0].size());

  Eigen::Vector2d center;
  center.x() = 0;
  center.y() = 0;

  G.resize(vertexNumber, dim);
  h.resize(vertexNumber);
  G.setOnes();
  h.setOnes();

  for(auto & p : inputPoints)
  {
    center.x() += p.x();
    center.y() += p.y();
  }

  center.x() = center.x() / (double)vertexNumber;
  center.y() = center.y() / (double)vertexNumber;

  int vNumber = 0;

  for(auto idx = inputPoints.begin(); idx != inputPoints.end(); idx++, vNumber++)
  {
    Eigen::Vector2d point_one;
    point_one.x() = idx->x();
    point_one.y() = idx->y();

    Eigen::Vector2d point_two;

    if((idx + 1) == inputPoints.end())
    {
      point_two.x() = inputPoints.begin()->x();
      point_two.y() = inputPoints.begin()->y();
    }
    else
    {
      point_two.x() = (idx + 1)->x();
      point_two.y() = (idx + 1)->y();
    }

    Eigen::Vector2d difference = point_two - point_one;
    // difference.normalize();
    double slope = difference.y() / difference.x();

    clampSlope(slope, miniSlope, maxSlope);

    G(vNumber, 0) = -slope;
    h(vNumber) = -slope * point_one.x() + point_one.y();

    int lineSign = 1;
    /// should remove slope?
    if(!((center.y() - slope * center.x()) <= h(vNumber)))
    {
      lineSign = -1;
    }
    // Correct the sign with the centeroid point
    G.block(vNumber, 0, 1, 2) *= lineSign;
    h(vNumber) *= lineSign;

  } // end of iterating over points

} // end of pointsToInequalityMatrix

template void mc_impact::pointsToInequalityMatrix<Eigen::Vector2d>(const std::vector<Eigen::Vector2d> & inputPoints,
                                                                   Eigen::MatrixXd & G,
                                                                   Eigen::VectorXd & h,
                                                                   double miniSlope,
                                                                   double maxSlope);

template void mc_impact::pointsToInequalityMatrix<Eigen::Vector3d>(const std::vector<Eigen::Vector3d> & inputPoints,
                                                                   Eigen::MatrixXd & G,
                                                                   Eigen::VectorXd & h,
                                                                   double miniSlope,
                                                                   double maxSlope);

template void mc_impact::pointsToInequalityMatrix<StaticPoint>(const std::vector<StaticPoint> & inputPoints,
                                                               Eigen::MatrixXd & G,
                                                               Eigen::VectorXd & h,
                                                               double miniSlope,
                                                               double maxSlope);

template<typename Point>
void mc_impact::pointsToInequalityMatrix(const std::vector<std::shared_ptr<Point>> & inputPoints,
                                         Eigen::MatrixXd & G,
                                         Eigen::VectorXd & h,
                                         double miniSlope,
                                         double maxSlope)
{
  int vertexNumber = static_cast<int>(inputPoints.size());
  int dim = static_cast<int>(inputPoints[0]->size());

  Eigen::Vector2d center;
  center.x() = 0;
  center.y() = 0;

  G.resize(vertexNumber, dim);
  h.resize(vertexNumber);
  G.setOnes();
  h.setOnes();

  for(auto & p : inputPoints)
  {
    center.x() += p->x();
    center.y() += p->y();
  }
  center.x() = center.x() / (double)vertexNumber;
  center.y() = center.y() / (double)vertexNumber;

  int vNumber = 0;

  for(auto idx = inputPoints.begin(); idx != inputPoints.end(); idx++, vNumber++)
  // for(auto & idx: inputPoints)
  {
    Eigen::Vector2d point_one;
    point_one.x() = (*idx)->x();
    point_one.y() = (*idx)->y();

    Eigen::Vector2d point_two;

    if((idx + 1) == inputPoints.end())
    {
      point_two.x() = (*inputPoints.begin())->x();
      point_two.y() = (*inputPoints.begin())->y();
    }
    else
    {
      point_two.x() = (*(idx + 1))->x();
      point_two.y() = (*(idx + 1))->y();
    }

    Eigen::Vector2d difference = point_two - point_one;
    // difference.normalize();
    double slope = difference.y() / difference.x();

    clampSlope(slope, miniSlope, maxSlope);

    G(vNumber, 0) = -slope;
    h(vNumber) = -slope * point_one.x() + point_one.y();

    int lineSign = 1;
    /// should remove slope?
    if(!((center.y() - slope * center.x()) <= h(vNumber)))
    {
      lineSign = -1;
    }
    // Correct the sign with the centeroid point
    G.block(vNumber, 0, 1, 2) *= lineSign;
    h(vNumber) *= lineSign;

  } // end of iterating over points
  vNumber++;
} // end of pointsToInequalityMatrix

template void mc_impact::pointsToInequalityMatrix<Eigen::Vector2d>(
    const std::vector<std::shared_ptr<Eigen::Vector2d>> & inputPoints,
    Eigen::MatrixXd & G,
    Eigen::VectorXd & h,
    double miniSlope,
    double maxSlope);

template void mc_impact::pointsToInequalityMatrix<Eigen::Vector3d>(
    const std::vector<std::shared_ptr<Eigen::Vector3d>> & inputPoints,
    Eigen::MatrixXd & G,
    Eigen::VectorXd & h,
    double miniSlope,
    double maxSlope);

template void mc_impact::pointsToInequalityMatrix<StaticPoint>(
    const std::vector<std::shared_ptr<StaticPoint>> & inputPoints,
    Eigen::MatrixXd & G,
    Eigen::VectorXd & h,
    double miniSlope,
    double maxSlope);
