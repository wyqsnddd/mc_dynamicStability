#pragma once

#include <mc_rbdyn/Robots.h>

#include "McDynamicStability/McContact.h"
#include "McDynamicStability/McPolytopeDescriptor.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <eigen-lssol/LSSOL_QP.h>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <polytope/staticStabilityPolytope.h>

namespace mc_impact
{

struct IeqConstraintBlocks
{
  Eigen::MatrixXd G_zmp;
  Eigen::VectorXd h_zmp;
};

struct McZMPAreaParams
{
  unsigned iterationLimit = 50;
  double convergeThreshold = 0.01;
};

template<typename Point>
class McZMPArea
/*! \brief This is an c++ implementation of the multi-contact-ZMP-area calculation
 *  1. We use  the "Ray-shooting-method" by `Bretl and Lall's algorithm`.
 *  2. Following the examples by: `project_polytope_bretl() in pymanoid`.
 *  3. Returns (*) verticies of the `multi-contact-zmp-area`.
 */

{

  static constexpr double LOWER_SLOPE = 0.01;
  static constexpr double UPPER_SLOPE = 100.0;

  ///< Size of the Grasp Matrix.
  static constexpr int GM_SIZE = 6;

  ///< Size of the Rotation Matrix.
  static constexpr int RM_SIZE = 3;

public:
  McZMPArea(const mc_rbdyn::Robot & robot,
            std::shared_ptr<McContactSet> contactSetPtr,
            const McZMPAreaParams & mcZMPAreaParams);
  ~McZMPArea() {}

  /*! It needs to be updated in each iteration.
   * \param  zmpVerticies saves the updated verticies of the ZMP
   * \param  height of the surface where the ZMP is projected. The default value is 0.0.
   */
  void computeMcZMPArea(double height = 0.0);

  /*! Obtain a reference to the robot.
   */
  inline const mc_rbdyn::Robot & getRobot()
  {
    return robot_;
  }

  ///< Get pointer to the set of contacts.
  inline std::shared_ptr<McContactSet> getContactSet()
  {
    return contactsPtr_;
  }

  inline const IeqConstraintBlocks & getIeqConstraint()
  {
    return ieqConstraintBlocks_;
  }

  inline int getNumVertex() const
  {
    return numVertex_;
  }

  inline int getMaxNumVertex() const
  {
    return polytopeProjectorPtr_->getMaxIteration();
  }
  inline const std::vector<Eigen::Vector2d> & getPolygonVertices()
  {
    return polygonVertices_;
  }
  const McZMPAreaParams & getParams() const
  {
    return McZMPAreaParams_;
  }

private:
  const mc_rbdyn::Robot & robot_;

  std::shared_ptr<McPolytopeDescriptor> pdPtr_; ///< The interface object for the stabiliPlus librarry

  int numVertex_ = 0;

  std::shared_ptr<StaticStabilityPolytope> polytopeProjectorPtr_;

  std::shared_ptr<McContactSet> contactsPtr_;
  void update_();

  McZMPAreaParams McZMPAreaParams_;

  IeqConstraintBlocks ieqConstraintBlocks_;

  std::vector<Eigen::Vector2d> polygonVertices_;

  // std::map<std::string, std::McContact> contacts_;
};

} // namespace mc_impact
