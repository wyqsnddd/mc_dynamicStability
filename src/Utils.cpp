#include "McDynamicStability/Utils.h"

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullPoints.h>

//#include <polytope/staticPoint.h>

Eigen::Matrix3d mc_impact::crossMatrix(const Eigen::Vector3d & input)
{

  Eigen::Matrix3d skewSymmetricMatrix = Eigen::Matrix3d::Zero();

  skewSymmetricMatrix(0, 1) = -input(2);
  skewSymmetricMatrix(1, 0) = input(2);

  skewSymmetricMatrix(0, 2) = input(1);
  skewSymmetricMatrix(2, 0) = -input(1);

  skewSymmetricMatrix(1, 2) = -input(0);
  skewSymmetricMatrix(2, 1) = input(0);

  return skewSymmetricMatrix;
}

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

/*
template void mc_impact::pointsToInequalityMatrix<StaticPoint>(const std::vector<StaticPoint> & inputPoints,
                                                               Eigen::MatrixXd & G,
                                                               Eigen::VectorXd & h,
                                                               double miniSlope,
                                                               double maxSlope);

                     */
template<typename Point>
void mc_impact::pointsToInequalityMatrix(const std::vector<Point> & inputPoints,
                                         Eigen::MatrixXd & G,
                                         Eigen::VectorXd & h,
                                         std::vector<Eigen::Vector2d> & pointsOut,
                                         double miniSlope,
                                         double maxSlope)
{

  int points_size = static_cast<int>(inputPoints.size());
  // std::cerr<< "The input points size is: " << points_size << std::endl;
  // pointsOut.resize(points_size);

  int dim = static_cast<int>(inputPoints[0].size());
  // Go through the points: update the verticies_
  std::vector<double> points_in;
  points_in.reserve(points_size * dim);

  Eigen::Vector2d center;
  center.x() = 0;
  center.y() = 0;

  for(const auto & p : inputPoints)
  {
    points_in.push_back(p.x());
    points_in.push_back(p.y());
    // pointsOut.emplace_back(p->x(), p->y());

    center.x() += p.x();
    center.y() += p.y();
  }

  center.x() = center.x() / (double)points_size;
  center.y() = center.y() / (double)points_size;

  // std::cerr "The center is: " << center.transpose() << std::endl;

  // Preprocess the points: order the vertices.
  orgQhull::Qhull qhull;
  // std::cout << "The Qhull dim is: " << dim << std::endl;

  qhull.runQhull("", dim, points_size, points_in.data(), "Qt");

  G.resize(points_size, dim);
  h.resize(points_size);
  G.setOnes();
  h.setOnes();

  // std::cerr<<"The center is "<<  qhull.feasiblePoint()[0]<<", "<< qhull.feasiblePoint()[1]<<std::endl;
  // std::cerr "The vertices of the porjected polygon is: " << std::endl;

  // std::cerr<< "------------------------------------------ " << std::endl;

  auto tempVertexList = qhull.vertexList();
  // std::cout<<"The center is "<<  qhull.feasiblePoint()[0]<<", "<< qhull.feasiblePoint()[1]<<std::endl;
  // center << qhull.feasiblePoint()[0], qhull.feasiblePoint()[1]; //, qhull.feasiblePoint()[2];

  // int qhullVertexNumer = qhull.vertexCount();
  // std::cout << "The qhull points num is: " << qhullVertexNumer << std::endl;
  // pointsOut.resize(qhullVertexNumer);

  int vNumber = 0;
  for(auto ii = tempVertexList.begin(); ii != tempVertexList.end(); ++ii, ++vNumber)
  {

    /*
    std::cerr<< "Processing qhull point " << ii->point().coordinates()[0] << ", " << ii->point().coordinates()[1]
              << std::endl;
        */

    // pointsOut.push_back({ii->point().coordinates()[0], ii->point().coordinates()[1]});
    pointsOut.emplace_back(ii->point().coordinates()[0], ii->point().coordinates()[1]);

    Eigen::Vector2d point_one, point_two;
    point_one << ii->point().coordinates()[0], ii->point().coordinates()[1];

    if((ii + 1) == tempVertexList.end())
    {
      auto jj = tempVertexList.begin();
      point_two << jj->point().coordinates()[0], jj->point().coordinates()[1];
    }
    else
    {
      point_two << (ii + 1)->point().coordinates()[0], (ii + 1)->point().coordinates()[1];
    }

    Eigen::Vector2d difference = point_two - point_one;
    // difference.normalize();
    double slope = difference.y() / difference.x();

    clampSlope(slope, miniSlope, maxSlope);

    G(vNumber, 0) = -slope;
    h(vNumber) = -slope * point_one.x() + point_one.y();

    int lineSign = 1;
    if(!(center.y() - slope * center.x() <= -slope * (h(vNumber))))
    {
      lineSign = -1;
    }
    // Correct the sign with the centeroid point
    G.block(vNumber, 0, 1, 2) *= lineSign;
    h(vNumber) *= lineSign;
  } // iterate for each vertex

} // end of pointsToInequalityMatrix

template void mc_impact::pointsToInequalityMatrix<Eigen::Vector2d>(const std::vector<Eigen::Vector2d> & inputPoints,
                                                                   Eigen::MatrixXd & G,
                                                                   Eigen::VectorXd & h,
                                                                   std::vector<Eigen::Vector2d> & points,
                                                                   double miniSlope,
                                                                   double maxSlope);

template void mc_impact::pointsToInequalityMatrix<Eigen::Vector3d>(const std::vector<Eigen::Vector3d> & inputPoints,
                                                                   Eigen::MatrixXd & G,
                                                                   Eigen::VectorXd & h,
                                                                   std::vector<Eigen::Vector2d> & points,
                                                                   double miniSlope,
                                                                   double maxSlope);
/*
template void mc_impact::pointsToInequalityMatrix<StaticPoint>(
    const std::vector<StaticPoint> & inputPoints,
    Eigen::MatrixXd & G,
    Eigen::VectorXd & h,
    std::vector<Eigen::Vector2d> & points,
    double miniSlope,
    double maxSlope);
    */

template<typename Point>
void mc_impact::pointsToInequalityMatrixSimple(const std::vector<Point> & inputPoints,
                                               Eigen::MatrixXd & G,
                                               Eigen::VectorXd & h,
                                               double miniSlope,
                                               double maxSlope)
{

  int points_size = static_cast<int>(inputPoints.size());
  // std::cerr<< "The input points size is: " << points_size << std::endl;
  // pointsOut.resize(points_size);

  int dim = static_cast<int>(inputPoints[0].size());
  // Go through the points: update the verticies_
  Eigen::Vector2d center;
  center.x() = 0;
  center.y() = 0;

  for(const auto & p : inputPoints)
  {
    center.x() += p.x();
    center.y() += p.y();
  }

  center.x() = center.x() / (double)points_size;
  center.y() = center.y() / (double)points_size;

  // std::cerr "The center is: " << center.transpose() << std::endl;

  G.resize(points_size, dim);
  h.resize(points_size);
  G.setOnes();
  h.setOnes();

  // std::cerr<<"The center is "<<  qhull.feasiblePoint()[0]<<", "<< qhull.feasiblePoint()[1]<<std::endl;
  // std::cerr "The vertices of the porjected polygon is: " << std::endl;

  // std::cerr<< "------------------------------------------ " << std::endl;

  // std::cout<<"The center is "<<  qhull.feasiblePoint()[0]<<", "<< qhull.feasiblePoint()[1]<<std::endl;
  // center << qhull.feasiblePoint()[0], qhull.feasiblePoint()[1]; //, qhull.feasiblePoint()[2];

  // int qhullVertexNumer = qhull.vertexCount();
  // std::cout << "The qhull points num is: " << qhullVertexNumer << std::endl;
  // pointsOut.resize(qhullVertexNumer);

  int vNumber = 0;
  for(auto ii = inputPoints.begin(); ii != inputPoints.end(); ++ii, ++vNumber)
  {

    /*
    std::cerr<< "Processing qhull point " << ii->point().coordinates()[0] << ", " << ii->point().coordinates()[1]
              << std::endl;
        */

    // pointsOut.push_back({ii->point().coordinates()[0], ii->point().coordinates()[1]});

    Eigen::Vector2d point_one, point_two;
    point_one << ii->x(), ii->y();

    if((ii + 1) == inputPoints.end())
    {
      auto jj = inputPoints.begin();
      point_two << jj->x(), jj->y();
    }
    else
    {
      point_two << (ii + 1)->x(), (ii + 1)->y();
    }

    Eigen::Vector2d difference = point_two - point_one;
    // difference.normalize();
    double slope = difference.y() / difference.x();

    clampSlope(slope, miniSlope, maxSlope);

    G(vNumber, 0) = -slope;
    h(vNumber) = -slope * point_one.x() + point_one.y();

    int lineSign = 1;
    if(!(center.y() - slope * center.x() <= -slope * (h(vNumber))))
    {
      lineSign = -1;
    }
    // Correct the sign with the centeroid point
    G.block(vNumber, 0, 1, 2) *= lineSign;
    h(vNumber) *= lineSign;
  } // iterate for each vertex

} // end of pointsToInequalityMatrix

template void mc_impact::pointsToInequalityMatrixSimple<Eigen::Vector2d>(
    const std::vector<Eigen::Vector2d> & inputPoints,
    Eigen::MatrixXd & G,
    Eigen::VectorXd & h,
    double miniSlope,
    double maxSlope);

template void mc_impact::pointsToInequalityMatrixSimple<Eigen::Vector3d>(
    const std::vector<Eigen::Vector3d> & inputPoints,
    Eigen::MatrixXd & G,
    Eigen::VectorXd & h,
    double miniSlope,
    double maxSlope);

void mc_impact::removeDuplicates(std::vector<Eigen::Vector2d> & vec)
{

  // If we would like to store values, we should use std::unordered_map.
  std::unordered_set<Eigen::Vector2d, mc_impact::ApproxHash, mc_impact::ApproxEqual> pointSet;

  auto ii = vec.begin();
  while(ii != vec.end())
  {

    // std::cerr<<yellow<<"Processing: "<<ii->transpose()<<reset<<std::endl;
    if(pointSet.find(*ii) != pointSet.end()) // O(1) lookup time for unordered_set
    {

      // std::cerr<<red<<"Found duplicate: "<<ii->transpose()<<reset<<std::endl;
      vec.erase(ii); // vec.erase returns the next valid iterator
    }
    else
    {
      pointSet.insert(*ii);
      // std::cerr<<green<<"Inserted: "<<ii->transpose()<<reset<<std::endl;

      // std::cout<<"pointset size: "<<pointset.size()<<std::endl;

      /*
      for(auto & jj : pointset){
       std::cerr << jj.transpose()<< std::endl;
      }
      */
      // std::cout<<"If inserted: "<<(pointset.find(*ii) == pointset.end())<<std::endl;
      ii++;
    }
  }

  // std::cout<<red<<"unordered pointset size: "<<pointSet.size()<<magenta<<std::endl;

} // end of removeDuplicates

/*

bool pointsAreClose(const Eigen::Vector2d & pointOne, const Eigen::Vector2d & pointTwo)
{

   double threshold = 0.01;
   if( (fabs(pointOne.x() - pointTwo.x())<threshold )&& (fabs(pointOne.y() - pointTwo.y()) < threshold))
     return true;
   else
     return false;
}
void mc_impact::removeDuplicatesSimple(std::vector<Eigen::Vector2d>& inputPoints)
{

  int vNumber = 0;
  for(auto ii = inputPoints.begin(); ii != inputPoints.end(); ++ii, ++vNumber)
  {
    Eigen::Vector2d point_one, point_two;
    point_one << ii->x(), ii->y();
    if((ii + 1) == inputPoints.end())
    {
      auto jj = inputPoints.begin();
      point_two << jj->x(), jj->y();
    }
    else
    {
      point_two << (ii + 1)->x(), (ii + 1)->y();
    }
    if(pointsAreClose(point_one, point_two)){
      std::cout<<red<<"removed "<<(ii+1)->transpose()<<reset<<std::endl;
      inputPoints.erase(ii+1);
    }

  }


} // end of removeDuplicates

*/

// template void mc_impact::removeDuplicates<Eigen::Vector2d>(std::vector<Eigen::Vector2d>& vec);
