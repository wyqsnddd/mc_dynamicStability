#pragma once

#include <mc_control/fsm/Controller.h>
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
  bool useSpatialVectorAlgebra = false; ///< Use the sva-consistent representation: i.e. wrench = [\tau, f], otherwise
                                        ///< we represent wrench = [f, \tau].
  bool initialContactStatus;
};

class McContact
{
public:
  McContact(const McContactParams & inputParams, const mc_rbdyn::Robot & robot);

  ~McContact() {}

  /*! \brief Updates: (1) CoP (2) Contact Vertices (3) Grasp Matrix.
   */
  //void update(const mc_rbdyn::Robot & robot);
  void update();

  inline const McContactParams & getContactParams() const
  {
    return mcContactParams_;
  }

  /*! \brief get the contact-wrench cone in the local contact frame.
   */
  inline const Eigen::MatrixXd & contactWrenchCone() const
  {
    return CWC_;
  }

  /*! \brief get the contact-wrench cone in the inertial contact frame.
   */
  inline const Eigen::MatrixXd & contactWrenchConeInertialFrame() const
  {
    return CWCInertial_; 
  }

  inline const sva::ForceVecd & desiredWrench()
  {
    return desiredWrench_;
  }

  /*! \brief desired CoP in the contact surface frame, e.g. LeftFoot or RightFoot.
   */
  inline const Eigen::Vector3d & desiredCop()
  {
    return desiredCoP_;
  }

  /*! \brief returns the vertices of the contact area in the INERTIAL frame.
   */
  const std::vector<Eigen::Vector3d> & getContactAreaVertices() const
  {
    return inertialContactAreaVertices_;
  }

  /*! \brief returns the vertices of the surface in the INERTIAL frame.
   * If the planar contact is not idea, we assume that the surface area is smaller than the contact area.
   */
  const std::vector<Eigen::Vector3d> & getSurfaceVertices() const
  {
    return inertialSurfaceVertices_;
  }

  /*! \brief MeasuredCoP of the contact surface frame, e.g. LeftFoot or RightFoot, in the INERTIAL frame.
   */
  inline const Eigen::Vector3d & measuredCop() const
  {
    return measuredCoP_;
  }

  /*! \brief Wrench*Multiplier gives the resultant wrench at the origin of the intertial frame
   */
  inline const Eigen::Matrix6d & getResultantWrenchMultiplier() const
  {
    return resultantWrenchMultiplier_; 
  }
  
  inline const Eigen::Matrix6d & getGraspMatrix() const
  {
    return graspMatrix_;
  }
  inline bool useSVA() const
  {
    return getContactParams().useSpatialVectorAlgebra;
  }

  inline void breakContact() 
  {
    inContact_ = false; 
  }
  inline void setContact() 
  {
    inContact_ = true; 
  }
  /*! \return if the contact is set 
   */
  inline bool inContact() const
  {
     return inContact_; 
  }

  /*! Add gui items (contact points, surface normals and eg. ) to rviz
   */
  void addGuiItems(mc_control::fsm::Controller &ctl) const;


  /*! \return reference to the robot (either pyhsices-engine simulated or real)
   */
  inline const mc_rbdyn::Robot & robot() const
  {
    return robot_;
  }

private:
  McContactParams mcContactParams_;

  const mc_rbdyn::Robot & robot_;

  sva::ForceVecd desiredWrench_ = sva::ForceVecd::Zero(); ///< Desired wrench in the body frame, e.g. l/r_sole ;
  Eigen::Vector3d desiredCoP_ =
      Eigen::Vector3d::Zero(); ///< Desired CoP in the contact surface frame, e.g. LeftFoot or RightFoot.
  Eigen::Vector3d measuredCoP_ =
      Eigen::Vector3d::Zero(); ///< MeasuredCoP in the contact surface frame, e.g. LeftFoot or RightFoot.
  Eigen::MatrixXd CWC_; ///< Contact wrench cone in the local contact frame, e.g. l/r_sole ;
  Eigen::MatrixXd CWCInertial_; ///< Contact wrench cone in the inertial contact frame, i.e. the robot inertial frame: X_0.
  Eigen::Matrix6d resultantWrenchMultiplier_; ///< Wrench*Multiplier gives the resultant wrench at the origin of the intertial frame

  void updateContactAreaVerticiesAndCoP_();
  //void updateContactAreaVerticiesAndCoP_(const mc_rbdyn::Robot & robot);

  //void updateCWCInertial_(const mc_rbdyn::Robot & robot);
  void updateCWCInertial_();

  inline const std::vector<Eigen::Vector3d> & getLocalContactAreaVerticies_() const
  {
    return localContactAreaVertices_;
  }
  std::vector<Eigen::Vector3d> localContactAreaVertices_;

  std::vector<Eigen::Vector3d> inertialContactAreaVertices_;
  std::vector<Eigen::Vector3d> inertialSurfaceVertices_;

  /*! \brief Update the measured CoP in each iteration according to the input measured Wrench.
   * \param input Wrench.
   */
  // void updateCoP_(const mc_rbdyn::Robot & robot);

  /*
   * \brief calculate the grasp matrix w.r.t. the inertial frame
   *  In the future, we need one more argument: the name of the target frame.
   *  this function follows the `spatial-vector-algebra`.
   */
  void calcGraspMatrix_();

  /* \brief This function follows the `GeometricRobotics`.
   * The difference is that we assume the wrench is given as: [Torque, Force].
   */
  void calcGeometricGraspMatrix_();
  //void calcGeometricGraspMatrix_(const mc_rbdyn::Robot & robot);

  void initializeCWC_();

  Eigen::Matrix6d graspMatrix_ = Eigen::Matrix6d::Identity();

  bool inContact_ = false;
};

struct McContactSet
{
public:
  McContactSet() {}
  ~McContactSet() {}
  /*! Obtain a reference to the contact.
   */
  const McContact & getContact(const std::string & name) const;

  bool addContact(const McContactParams & inputParams, const mc_rbdyn::Robot & robot);

  inline const std::map<std::string, mc_impact::McContact> & getContactMap() const
  {
    return contacts_;
  }

  //void update(const mc_rbdyn::Robot & robot);
  void update();

  /*! Add gui items (contact points, surface normals and eg. ) of all the contacts to rviz
   */
  void addGuiItems(mc_control::fsm::Controller &ctl) const;

private:
  std::map<std::string, mc_impact::McContact> contacts_;
};

} // namespace mc_impact
