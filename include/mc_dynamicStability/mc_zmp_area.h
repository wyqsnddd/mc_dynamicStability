#pragma once

#include <mc_rbdyn/Robots.h>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <eigen-lssol/LSSOL_QP.h>
#include <iomanip>
#include <iostream>
#include <math.h>

#include <polytope/staticStabilityPolytope.h>

#include "mc_dynamicStability/mc_polytopeDescriptor.h"
#include "mc_dynamicStability/contact.h"

namespace mc_impact { 
struct McZMPAreaParams
{

  double x;
};

template<typename Point>
class McZMPArea
/*! \brief This is an c++ implementation of the multi-contact-ZMP-area calculation
 *  1. We use  the "Ray-shooting-method" by `Bretl and Lall's algorithm`.
 *  2. Following the examples by: `project_polytope_bretl() in pymanoid`.
 *  3. Returns (*) verticies of the `multi-contact-zmp-area`.
 */

{

public:
  McZMPArea(const mc_rbdyn::Robot & robot, const struct McZMPAreaParams params);
  ~McZMPArea() {}

  /*! It needs to be updated in each iteration.
   * \param  zmpVerticies saves the updated verticies of the ZMP
   * \param  height of the surface where the ZMP is projected. The default value is 0.0.
   */
  void computeMcZMPArea(double height = 0.8);

  /*! Obtain a reference to the robot.
   */
  inline const mc_rbdyn::Robot & getRobot()
  {
    return robot_;
  }

  McContactSet contacts; ///< The set of contacts.

private:
  const mc_rbdyn::Robot & robot_;

  McZMPAreaParams params_;

  std::shared_ptr<McPolytopeDescriptor> pdPtr_; ///< The interface object for the stabiliPlus librarry
  
  std::shared_ptr<StaticStabilityPolytope> polytopeProjectorPtr_;

  void update_();

  Eigen::Matrix3d crossMatrix_(const Eigen::Vector3d & input);

  // std::map<std::string, std::McContact> contacts_;
};

} // namespace mc_impact
