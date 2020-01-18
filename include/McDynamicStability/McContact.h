#pragma once

#include <mc_rbdyn/Robots.h>

#include "Utils.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace mc_impact
{

struct McContactParams
/*!
 * \brief Parameters of a contact. We use the surface name to denote this contact
 */
{
  std::string surfaceName; ///< surface name of the contact.
  std::string bodyName; ///< Body name of the link where the contact is defined.
  std::string sensorName; ///< Sensor name of the link.

  double halfX = 0.06;
  double halfY = 0.05;
  double frictionCoe = 0.7; ///< Friction coefficient.
  double minForce = 15.0; ///< The minimum contact force.
  // double maxPressure = 1000.0; ///< The minimum contact force.
  int index = 0; ///< Index of the contact in the wrenchDistributionQP.
  bool useSpatialVectorAlgebra = false;///< Use the sva-consistent representation: i.e. wrench = [\tau, f], otherwise we represent wrench = [f, \tau].
};

class McContact
{
public:
  McContact(const McContactParams & inputParams, const mc_rbdyn::Robot & robot);

  ~McContact() {}

  /*! \brief Updates: (1) CoP (2) Contact Vertices (3) Grasp Matrix.
   */
  void update(const mc_rbdyn::Robot & robot);

  inline const McContactParams & getContactParams() const
  {
    return mcContactParams_;
  }

  inline const Eigen::MatrixXd & contactWrenchCone() const
  {
    return CWC_;
  }

  inline const sva::ForceVecd & desiredWrench()
  {
    return desiredWrench_;
  }

  /*! \brief desired CoP in the contact surface frame, e.g. LeftFoot or RightFoot.
  */
  inline const Eigen::Vector2d & desiredCop()
  {
    return desiredCoP_;
  }

  /*! \brief returns the vertices in the inertial frame.
   */
  const std::vector<Eigen::Vector3d> & getContactAreaVerticies()
  {
    return inertialContactAreaVertices_;
  }

  /*! \brief MeasuredCoP in the contact surface frame, e.g. LeftFoot or RightFoot.
  */
  inline const Eigen::Vector2d & measuredCop()
  {
    return measuredCoP_;
  }

  inline const Eigen::Matrix6d & getGraspMatrix() const
  {
    return graspMatrix_; 
  }
  inline bool useSVA() const
  {
   return getContactParams().useSpatialVectorAlgebra; 
  }
private:
  McContactParams mcContactParams_;

  sva::ForceVecd desiredWrench_= sva::ForceVecd::Zero(); ///< Desired wrench in the body frame, e.g. l/r_sole ;
  Eigen::Vector2d desiredCoP_ = Eigen::Vector2d::Zero(); ///< Desired CoP in the contact surface frame, e.g. LeftFoot or RightFoot.
  Eigen::Vector2d measuredCoP_= Eigen::Vector2d::Zero(); ///< MeasuredCoP in the contact surface frame, e.g. LeftFoot or RightFoot.
  Eigen::MatrixXd CWC_; ///< Contact wrench cone in the local contact frame, e.g. l/r_sole ;


  void updateContactAreaVerticies_(const mc_rbdyn::Robot & robot);

  inline const std::vector<Eigen::Vector3d> & getLocalContactAreaVerticies_() const
  {
    return localContactAreaVertices_; 
  }
  std::vector<Eigen::Vector3d> localContactAreaVertices_;

  std::vector<Eigen::Vector3d> inertialContactAreaVertices_;


  /*! \brief Update the measured CoP in each iteration according to the input measured Wrench. 
   * \param input Wrench. 
  */
  void updateCoP_(const sva::ForceVecd & inputWrench);

    /*
   * \brief calculate the grasp matrix w.r.t. the inertial frame
   *  In the future, we need one more argument: the name of the target frame.
   *  this function follows the `spatial-vector-algebra`.
   */
  void calcGraspMatrix_(const mc_rbdyn::Robot & robot);

  /* \brief This function follows the `GeometricRobotics`.
   * The difference is that we assume the wrench is given as: [Torque, Force].
   */
  void calcGeometricGraspMatrix_(const mc_rbdyn::Robot & robot);



  void initializeCWC_();


  Eigen::Matrix6d graspMatrix_ = Eigen::Matrix6d::Identity();
};

struct McContactSet
{
public:
  McContactSet() {}
  ~McContactSet() {}
  /*! Obtain a reference to the contact.
   */
  const McContact & getContact(const std::string & name);

  bool addContact(const McContactParams & inputParams, const mc_rbdyn::Robot & robot);

  inline const std::map<std::string, mc_impact::McContact> & getContactMap() const
  {
    return contacts_;
  }

  void update(const mc_rbdyn::Robot & robot);

private:
  std::map<std::string, mc_impact::McContact> contacts_;
};

} // namespace mc_impact
