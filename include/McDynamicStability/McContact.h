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
};

class McContact
{
public:
  McContact(const McContactParams & inputParams);

  ~McContact() {}

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

  inline const Eigen::Vector2d & cop()
  {
    return desiredCoP_;
  }

  void updateCWC();
  void updateCoP(const sva::ForceVecd & inputWrench);

  /*
   * \brief calculate the grasp matrix w.r.t. the inertial frame
   *  In the future, we need one more argument: the name of the target frame.
   *  this function follows the `spatial-vector-algebra`.
   */
  void calcGraspMatrix(Eigen::Matrix6d & G, const mc_rbdyn::Robot & realRobot) const;

  /* \brief This function follows the `GeometricRobotics`.
   * The difference is that we assume the wrench is given as: [Torque, Force].
   */
  void calcGeometricGraspMatrix(Eigen::Matrix6d & G, const mc_rbdyn::Robot & realRobot) const;

private:
  McContactParams mcContactParams_;

  sva::ForceVecd desiredWrench_; ///< Desired wrench in the body frame, e.g. l/r_sole ;
  Eigen::Vector2d desiredCoP_; ///< Desired CoP in the contact surface frame, e.g. LeftFoot or RightFoot.
  Eigen::MatrixXd CWC_; ///< Contact wrench cone in the local contact frame, e.g. l/r_sole ;
};

struct McContactSet
{
public:
  McContactSet() {}
  ~McContactSet() {}
  /*! Obtain a reference to the contact.
   */
  const McContact & getContact(const std::string & name);

  bool addContact(const McContactParams & inputParams);

  inline const std::map<std::string, mc_impact::McContact> & getContactMap() const
  {
    return contacts_;
  }

private:
  std::map<std::string, mc_impact::McContact> contacts_;
};

} // namespace mc_impact
