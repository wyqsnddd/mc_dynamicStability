#pragma once

#include <mc_control/fsm/Controller.h>
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

#include <geos/version.h>

#include <geos/geom/CoordinateSequenceFactory.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/Polygon.h>

#include <geos/geom/Point.h>

namespace mc_impact
{

	/*
struct McComAreaParams
{
  unsigned iterationLimit = 50;
  double projectionRadius = 10.0;
  double convergeThreshold = 0.01;
  bool useLIPMAssumptions = true;
  bool debug = false;
};

*/
class McComArea
/*! \brief This is an c++ implementation to calculate the 'multi-contact Com static equilibrium area'
 *  1. We use  the "Ray-shooting-method" by `Bretl and Lall's algorithm`.
 *  2. NOT yet following any example. We ony use the equations on the paper.
 *  3. Returns (*) verticies of the `multi-contact-com-area`.
 */
{
  static constexpr double LOWER_SLOPE = 0.01;
  static constexpr double UPPER_SLOPE = 100.0;

  ///< Size of the Grasp Matrix.
  static constexpr int GM_SIZE = 6;

  ///< Size of the Rotation Matrix.
  static constexpr int RM_SIZE = 3;

public:
  McComArea(const mc_rbdyn::Robot & robot,
            std::shared_ptr<McContactSet> contactSetPtr,
            const McProjectionParams & mcProjectionParams);
  ~McComArea() {}

  /*! It needs to be updated in each iteration.
   * \param  height of the surface where the ZMP is projected. The default value is 0.0.  */
  // void computeMcZMPArea(double height);
  void updateMcComArea();

  /*! Obtain a reference to the robot.
   */
  inline const mc_rbdyn::Robot & getRobot() const
  {
    return robot_;
  }
  ///< Get pointer to the set of contacts.
  inline std::shared_ptr<McContactSet> getContactSet() const
  {
    return contactsPtr_;
  }

  /*
   inline const IeqConstraintBlocks & getIeqConstraint() const
   {
     return ieqConstraintBlocks_;
   }
   */
  inline int getNumVertex() const
  {
    return static_cast<int>(polygonVertices_.size());
  }
  inline int getMaxNumVertex() const
  {
    // Each iteration should generate a new vetex.
    //return polytopeProjectorPtr_->getMaxIteration();
    return static_cast<int>(getParams().iterationLimit);
  }

  /*! \brief Return vertices of the Com static equilibrium polytope
   *
   */
  inline const std::vector<Eigen::Vector2d> & getPolygonVertices() const
  {
    return polygonVertices_;
  }
  const McProjectionParams & getParams() const
  {
    return mcProjectionParams_;
  }

  inline const std::shared_ptr<StaticStabilityPolytope> getProjector() const
  {
    return polytopeProjectorPtr_;
  }

  /*!
   * \brief Get the centroid of the projected multi-contact com area. 
   */
  inline const Eigen::Vector3d getCentroid()const
  {
    return centroid_; 
  }

  void print() const
  {
    std::cerr<<red<<"McComArea has: "<< getNumVertex()<<" vertices (maximum: "<<getMaxNumVertex()<<reset<<std::endl; 
  }

  /*! Add gui items (Mc Com static equilibrium area) to rviz
   */
  void addGuiItems(mc_control::fsm::Controller &ctl) const;


private:
  const mc_rbdyn::Robot & robot_;

  std::shared_ptr<McPolytopeDescriptor> pdPtr_; ///< The interface object for the stabiliPlus librarry
  std::shared_ptr<StaticStabilityPolytope> polytopeProjectorPtr_;

  std::shared_ptr<McContactSet> contactsPtr_;

  void update_();
  void updateLIPMAssumptions_(int numContact, const Eigen::MatrixXd & inputG);
  void computeMcComArea_();

  McProjectionParams mcProjectionParams_;

  // IeqConstraintBlocks ieqConstraintBlocks_;

  std::vector<Eigen::Vector2d> polygonVertices_;

  void calcCentroid_();
  Eigen::Vector3d centroid_;
};
} // namespace mc_impact
