#include "McDynamicStability/McContact.h"

namespace mc_impact
{

McContact::McContact(const McContactParams & inputParams) : mcContactParams_(inputParams)
{

  updateCWC();
}

void McContact::calcGeometricGraspMatrix(Eigen::Matrix6d & G, const mc_rbdyn::Robot & realRobot) const
{

  G.setIdentity();
  auto X_0_c = realRobot.surfacePose(getContactParams().surfaceName);

  G.block<3, 3>(0, 0) = X_0_c.rotation();
  G.block<3, 3>(3, 3) = G.block<3, 3>(0, 0);

  G.block<3, 3>(0, 3) = crossMatrix(X_0_c.translation()) * X_0_c.rotation();
}
void McContact::calcGraspMatrix(Eigen::Matrix6d & G, const mc_rbdyn::Robot & realRobot) const
{
  sva::PTransformd X_c_0 = realRobot.surfacePose(getContactParams().surfaceName).inv();

  G = X_c_0.dualMatrix();
}

void McContact::updateCWC()
{
  double X = 2 * getContactParams().halfX;
  double Y = 2 * getContactParams().halfY;
  double mu = getContactParams().frictionCoe;

  CWC_.resize(16, 6);
  CWC_.setZero();

  // clang-format off
  CWC_ <<
      // mx,  my,  mz,  fx,  fy,            fz,
          0,   0,   0, -1,    0,           -mu, 
	  0,   0,   0, +1,    0,           -mu, 
	  0,   0,   0,  0,   -1,           -mu, 
	  0,   0,   0,  0,   +1,           -mu, 
	 -1,   0,   0,  0,    0,            -Y, 
	 +1,   0,   0,  0,    0,            -Y, 
	  0,  -1,   0,  0,    0,            -X, 
	  0,  +1,   0,  0,    0,            -X, 
	+mu, +mu,  -1, -Y,   -X, -(X + Y) * mu, 
	+mu, -mu,  -1, -Y,   +X, -(X + Y) * mu, 
	-mu, +mu,  -1, +Y,   -X, -(X + Y) * mu, 
	-mu, -mu,  -1, +Y,   +X, -(X + Y) * mu, 
	+mu, +mu,  +1, +Y,   +X, -(X + Y) * mu, 
	+mu, -mu,  +1, +Y,   -X, -(X + Y) * mu, 
	-mu, +mu,  +1, -Y,   +X, -(X + Y) * mu, 
	-mu, -mu,  +1, -Y,   -X, -(X + Y) * mu;
  // clang-format on 
}

const McContact & McContactSet::getContact(const std::string & name)
{
  auto opt = contacts_.find(name);
  if(opt != (contacts_.end()))
  {
    return opt->second;
  }
  else
  {
    throw std::runtime_error(std::string("McContactSet::getContact: '") + name + std::string("' is not found."));
  }
}

bool McContactSet::addContact(const McContactParams & inputParams)
{
  auto opt = contacts_.find(inputParams.surfaceName);


  if(opt != (contacts_.end()))
  {
    std::cout << "Endeffector " << inputParams.surfaceName << " already exists." << std::endl;
    return false;
  }


  contacts_.insert({inputParams.surfaceName, {inputParams}});

  //std::cout << "McContactSet: Adding end-effector contact: " << inputParams.surfaceName << std::endl;

  return true;
}

} // namespace mc_impact
