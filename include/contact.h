#pragma once

#include <mc_rbdyn/Robots.h>

#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace mc_impact
{

struct ContactParams
{
  std::string surfaceName; ///< surface name of the contact
  std::string bodyName;
  std::string sensorName;

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
  McContact(const ContactParams & inputParams);

  ~McContact() {}

  inline const ContactParams & getContactParams() const
  {
    return contactParams_;
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

private:
  ContactParams contactParams_;

  sva::ForceVecd desiredWrench_; ///< Desired wrench in the body frame, e.g. l/r_sole ;
  Eigen::Vector2d desiredCoP_; ///< Desired CoP in the contact surface frame, e.g. LeftFoot or RightFoot.
  Eigen::MatrixXd CWC_; ///< Contact wrench cone in the local contact frame, e.g. l/r_sole ;
};

struct McContactSet
{

  /*! Obtain a reference to the contact.
   */
  const McContact & getContact(const std::string & name);

  bool addContact(const ContactParams & inputParams);

  const std::map<std::string, McContact> & getContactMap();

private:
  std::map<std::string, McContact> contacts_;
};

} // namespace mc_impact
