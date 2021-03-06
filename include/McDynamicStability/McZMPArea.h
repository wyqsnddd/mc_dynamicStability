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

/*
struct McZMPAreaParams
{
  unsigned iterationLimit = 50;
  double projectionRadius = 10.0;
  double convergeThreshold = 0.01;
  bool useLIPMAssumptions = true;
  bool debug = false;
  bool useSpatialVectorAlgebra = false; ///< Use the sva-consistent representation: i.e. wrench = [\tau, f], otherwise
};
*/
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
            const McProjectionParams & mcProjectionParams);
  ~McZMPArea() {}

  /*! It needs to be updated in each iteration.
   * \param  height of the surface where the ZMP is projected. The default value is 0.0.  */
  // void computeMcZMPArea(double height);
  void updateMcZMPArea(double height);

  /*! Obtain a reference to the robot.
   */
  inline const mc_rbdyn::Robot & getRobot() const
  {
    return robot_;
  }

  ///< Get pointer to the set of contacts.
  inline std::shared_ptr<McContactSet> getContactSet() const
  {
    assert( contactsPtr_ != nullptr );

    return contactsPtr_;
  }

  inline const IeqConstraintBlocks & getIeqConstraint() const
  {
    return ieqConstraintBlocks_;
  }

  inline int getNumVertex() const
  {
    return static_cast<int>(polygonVertices_.size());
  }

  inline int getMaxNumVertex() const
  {
    // Each iteration should generate a new vetex.
    return getParams().iterationLimit;
    //return polytopeProjectorPtr_->getMaxIteration();
  }
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
    if(polytopeProjectorPtr_ == nullptr){
       throw std::runtime_error("PolytopeProjector is asked without being allocated yet. ");
    }else{
       return polytopeProjectorPtr_;
    }
  }

  /*! \brief Obtain the ZMP assuming multiple contacts.
   */
  const Eigen::Vector3d & getMcZMP() const
  {
    return mcZMP_;   
  }

  /*! \brief Obtain the ZMP assuming double feet support. 
   */
  
  const Eigen::Vector3d & getBipedalZMP() const
  {
    return bipedalZMP_;   
  }
  

  static Eigen::Vector3d zmpCalculation(const Eigen::Vector3d & normal, const sva::ForceVecd Wrench);

void print() const
  {
    std::cerr<<red<<"McZMPArea has: "<< getNumVertex()<<" vertices (maximum: "<<getMaxNumVertex()<<reset<<std::endl; 
  }

private:
  const mc_rbdyn::Robot & robot_;

  std::shared_ptr<McPolytopeDescriptor> pdPtr_; ///< The interface object for the stabiliPlus librarry

  int numVertex_ = 0;

  std::shared_ptr<StaticStabilityPolytope> polytopeProjectorPtr_;

  std::shared_ptr<McContactSet> contactsPtr_;
  void update_();

  void updateLIPMAssumptions_(int numContact, const Eigen::MatrixXd & inputG);

  void computeMcZMPArea_(double height);

  McProjectionParams mcProjectionParams_;

  Eigen::Vector3d mcZMP_;

  Eigen::Vector3d bipedalZMP_;

  void updateMcZMP_(); 
  void updateBipedalZMP_(); 

  IeqConstraintBlocks ieqConstraintBlocks_;

  std::vector<Eigen::Vector2d> polygonVertices_;

  // std::map<std::string, std::McContact> contacts_;
};

} // namespace mc_impact
