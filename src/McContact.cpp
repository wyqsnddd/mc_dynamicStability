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

  G.block<3, 3>(0, 0) = X_0_c.rotation().transpose();
  G.block<3, 3>(3, 3) = G.block<3, 3>(0, 0);

  G.block<3, 3>(3, 0) = -G.block<3, 3>(0, 0) * crossMatrix(X_0_c.translation());
  /*
  G<<
 -0.16214695, -0.98665114, -0.0150962 ,  0.        ,  0.        ,  0.        ,
  0.89721345, -0.14104597, -0.41846632,  0.        ,  0.        ,  0.        ,
  0.41075101, -0.08139755,  0.90810685,  0.        ,  0.        ,  0.        ,
  0.        ,  0.10849315, -0.00440193, -0.16214695, -0.98665114, -0.0150962 ,
  0.09865849,  0.        , -0.01626989,  0.89721345, -0.14104597, -0.41846632,
 -0.11977171,  0.00316472, -0.        ,  0.41075101, -0.08139755,  0.90810685;


   std::cout<<"Contact: "<<getContactParams().surfaceName<<", link number:"
     <<realRobot.bodyIndexByName(getContactParams().bodyName)<<std::endl; std::cout<<"Translation: "<<
     X_0_c.translation().transpose()<<std::endl; std::cout<<"Rotation: "<<std::endl<< X_0_c.rotation()<<std::endl;

   Eigen::Quaterniond q(X_0_c.rotation());
   std::cout<<"Quaternion: "<<q.coeffs()<<std::endl;

   //Eigen::Vector3d euler = X_0_c.rotation().eulerAngles(2, 1, 0);

   //std::cout<<"yaw: "<< euler(0, 0)<<", pitch: "<<euler(1,0)<<", roll: "<<euler(3,0)<<std::endl;

   std::cout<<"The robot com position: "<< realRobot.com().transpose()<<", mass"<<realRobot.mass()<<std::endl;
*/
}

void McContact::calcGraspMatrix(Eigen::Matrix6d & G, const mc_rbdyn::Robot & realRobot) const
{
  sva::PTransformd X_c_0 = realRobot.surfacePose(getContactParams().surfaceName).inv();

  G = X_c_0.dualMatrix();

  G.block<3, 3>(3, 0) = G.block<3, 3>(0, 3);

  G.block<3, 3>(0, 3).setZero();
  // std::cout<<"The grasp matrix is: "<<std::endl<<G<<std::endl;
}

void McContact::updateCWC()
{
  double X = getContactParams().halfX;
  double Y = getContactParams().halfY;

  // Inner approximation
  double mu = getContactParams().frictionCoe / sqrt(2.0);

  CWC_.resize(16, 6);
  CWC_.setZero();

  // clang-format off
 // CWC_ <<
      // mx,  my,  mz,  fx,  fy,            fz,
      /*
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
	*/
  // clang-format on 
  
  // clang-format off
  CWC_ <<
      //  fx,  fy,            fz,  mx,  my,  mz,  
         -1,    0,           -mu,   0,   0,   0, 
	 +1,    0,           -mu,   0,   0,   0, 
	  0,   -1,           -mu,   0,   0,   0, 
	  0,   +1,           -mu,   0,   0,   0, 
	  0,    0,            -Y,  -1,   0,   0, 
	  0,    0,            -Y,  +1,   0,   0, 
	  0,    0,            -X,   0,  -1,   0, 
	  0,    0,            -X,   0,  +1,   0, 
	 -Y,   -X, -(X + Y) * mu, +mu, +mu,  -1, 
	 -Y,   +X, -(X + Y) * mu, +mu, -mu,  -1, 
	 +Y,   -X, -(X + Y) * mu, -mu, +mu,  -1, 
	 +Y,   +X, -(X + Y) * mu, -mu, -mu,  -1, 
	 +Y,   +X, -(X + Y) * mu, +mu, +mu,  +1, 
	 +Y,   -X, -(X + Y) * mu, +mu, -mu,  +1, 
	 -Y,   +X, -(X + Y) * mu, -mu, +mu,  +1, 
	 -Y,   -X, -(X + Y) * mu, -mu, -mu,  +1;
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
