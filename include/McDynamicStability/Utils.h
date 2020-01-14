#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <libqhullcpp/Qhull.h>
#include <memory>
//#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullPoints.h>
//#include <libqhullcpp/QhullVertexSet.h>

# include <unordered_set>
# include <cmath>
//# include <boost/container_hash/hash.hpp>

namespace mc_impact
{

const std::string red("\033[0;31m");
const std::string green("\033[1;32m");
const std::string yellow("\033[1;33m");
const std::string cyan("\033[0;36m");
const std::string magenta("\033[0;35m");
const std::string reset("\033[0m");

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


struct ApproxHash 
{
 std::size_t operator() (Eigen::Vector2d const& pt) const
 {


   size_t score = (size_t)(pt.x()*100) + (size_t)(pt.y()*10);
   std::cerr <<"Point: "<< pt.transpose()<< " has score: "<<score<<std::endl;
   return score; 
 }

};

struct ApproxEqual{
	// This is used to guarantee that no duplicates should happen when the hash collision happens. 
public:
 bool operator()(const Eigen::Vector2d & pt1, const Eigen::Vector2d & pt2) const {
    double threshold = 0.01;
    bool result = (fabs(pt1.x() - pt2.x())<threshold) && (fabs(pt1.y() - pt2.y())<threshold);

    //std::cerr<<cyan<<"Equal is called for: "<< pt1.transpose()<<" and "<<pt2.transpose()<<" which are " << result<<" equal. "<<" xdiff"<< xdiff<<", ydiff"<<ydiff<<reset<<std::endl;
    return result; 
}
};



//template<typename Point>
void removeDuplicates(std::vector<Eigen::Vector2d>& vec);
void removeDuplicatesNew(std::vector<Eigen::Vector2d>& vec);
void removeDuplicatesSimple(std::vector<Eigen::Vector2d>& vec);

template<typename Point>
void pointsToInequalityMatrix(const std::vector<Point> & inputPoints,
                              Eigen::MatrixXd & G,
                              Eigen::VectorXd & h,
                              double miniSlope = 0.01,
                              double maxSlope = 1000.0);

template<typename Point>
void pointsToInequalityMatrix(const std::vector<Point> & inputPoints,
                              Eigen::MatrixXd & G,
                              Eigen::VectorXd & h,
                              std::vector<Eigen::Vector2d> & point,
                              double miniSlope = 0.01,
                              double maxSlope = 1000.0);

template<typename Point>
void pointsToInequalityMatrixSimple(const std::vector<Point> & inputPoints,
                              Eigen::MatrixXd & G,
                              Eigen::VectorXd & h,
                              double miniSlope = 0.01,
                              double maxSlope = 1000.0);


Eigen::Matrix3d crossMatrix(const Eigen::Vector3d & input);

} // namespace mc_impact
