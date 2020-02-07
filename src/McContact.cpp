#include "McDynamicStability/McContact.h"

namespace mc_impact
{

McContact::McContact(const McContactParams & inputParams, const mc_rbdyn::Robot & robot) : mcContactParams_(inputParams)
{
  // Initialize the local contact points with a square:

  localContactAreaVertices_.emplace_back(getContactParams().halfX, getContactParams().halfY, 0.0);
  localContactAreaVertices_.emplace_back(getContactParams().halfX, -getContactParams().halfY, 0.0);
  localContactAreaVertices_.emplace_back(-getContactParams().halfX, -getContactParams().halfY, 0.0);
  localContactAreaVertices_.emplace_back(-getContactParams().halfX, getContactParams().halfY, 0.0);

  initializeCWC_();

  resultantWrenchMultiplier_.setIdentity();

#ifdef DEBUG
  std::cout << cyan << "Creating contact: " << getContactParams().surfaceName << reset << std::endl;
#endif

  update(robot);
}

void McContact::calcGeometricGraspMatrix_(const mc_rbdyn::Robot & robot)
{

  graspMatrix_.setIdentity();
  auto X_0_c = robot.surfacePose(getContactParams().surfaceName);

  // graspMatrix_.block<3, 3>(0, 0) = X_0_c.rotation().transpose();
  graspMatrix_.block<3, 3>(0, 0) = X_0_c.rotation();
  graspMatrix_.block<3, 3>(3, 3) = graspMatrix_.block<3, 3>(0, 0);

  // graspMatrix_.block<3, 3>(3, 0) = -graspMatrix_.block<3, 3>(0, 0) * crossMatrix(X_0_c.translation());

  auto rotation = X_0_c.rotation().transpose();
  // auto translation =  - rotation*X_0_c.translation();
  auto translation = X_0_c.translation();

  graspMatrix_.block<3, 3>(3, 0) = -rotation.transpose() * crossMatrix(translation);

  /*
  G<<
 -0.16214695, -0.98665114, -0.0150962 ,  0.        ,  0.        ,  0.        ,
  0.89721345, -0.14104597, -0.41846632,  0.        ,  0.        ,  0.        ,
  0.41075101, -0.08139755,  0.90810685,  0.        ,  0.        ,  0.        ,
  0.        ,  0.10849315, -0.00440193, -0.16214695, -0.98665114, -0.0150962 ,
  0.09865849,  0.        , -0.01626989,  0.89721345, -0.14104597, -0.41846632,
 -0.11977171,  0.00316472, -0.        ,  0.41075101, -0.08139755,  0.90810685;


   std::cout<<"Contact: "<<getContactParams().surfaceName<<", link number:"
     <<robot.bodyIndexByName(getContactParams().bodyName)<<std::endl; std::cout<<"Translation: "<<
     X_0_c.translation().transpose()<<std::endl; std::cout<<"Rotation: "<<std::endl<< X_0_c.rotation()<<std::endl;

   Eigen::Quaterniond q(X_0_c.rotation());
   std::cout<<"Quaternion: "<<q.coeffs()<<std::endl;

   //Eigen::Vector3d euler = X_0_c.rotation().eulerAngles(2, 1, 0);

   //std::cout<<"yaw: "<< euler(0, 0)<<", pitch: "<<euler(1,0)<<", roll: "<<euler(3,0)<<std::endl;

   std::cout<<"The robot com position: "<< robot.com().transpose()<<", mass"<<robot.mass()<<std::endl;
*/
}

void McContact::calcGraspMatrix_(const mc_rbdyn::Robot & robot)
{

  // --------------- Use the dual-matrix as the
  // sva::PTransformd X_c_0 = robot.surfacePose(getContactParams().surfaceName).inv();
  // graspMatrix_ = X_c_0.dualMatrix();

  // ---------------
  sva::PTransformd X_0_c = robot.surfacePose(getContactParams().surfaceName);
  graspMatrix_.setIdentity();
  /* This is a manual calculation */
  graspMatrix_.block<3, 3>(0, 0) = X_0_c.rotation();
  graspMatrix_.block<3, 3>(3, 3) = graspMatrix_.block<3, 3>(0, 0);

  auto rotation = X_0_c.rotation().transpose();
  // auto translation =  - rotation*X_0_c.translation();
  auto translation = X_0_c.translation();

  graspMatrix_.block<3, 3>(0, 3) = -rotation.transpose() * crossMatrix(translation);

  // std::cout<<"The grasp matrix is: "<<std::endl<<G<<std::endl;
}

void McContact::initializeCWC_()
{
  double X = 2 * getContactParams().halfX;
  double Y = 2 * getContactParams().halfY;

  // Inner approximation
  double mu = getContactParams().frictionCoe / sqrt(2.0);

  CWC_.resize(16, 6);
  CWC_.setZero();

  CWCInertial_.resize(16, 6);
  CWCInertial_.setZero();

  // clang-format off
  if(getContactParams().useSpatialVectorAlgebra)
  {
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
  }
  else
  {
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
  }
  // clang-format on 
  
  CWCInertial_ = CWC_;

}

void McContact::updateContactAreaVerticiesAndCoP_(const mc_rbdyn::Robot & robot)
{

  inertialContactAreaVertices_.clear();

  
  const sva::PTransformd & X_0_s = robot.surfacePose(getContactParams().surfaceName);
  //const sva::PTransformd & X_s_0 = robot.surfacePose(getContactParams().surfaceName).inv();

  // (0) Update the contact area points: 
  const auto & surface = robot.surface(getContactParams().surfaceName);
  const auto & pts = surface.points();
  auto rotation = X_0_s.rotation().inverse();
  //auto translation = - rotation * X_0_s.translation();
  auto translation = X_0_s.translation();




  std::cerr<<red<<"Updating contact: "<<getContactParams().surfaceName<<reset<<std::endl;
  std::cerr<<cyan<<"The surface: "<<surface.name()<<" has bodyName: "<<surface.bodyName()<<reset<<std::endl;

  std::cerr<<yellow<<"Rotation: "<<std::endl<<rotation<<reset<<std::endl;
  std::cerr<<yellow<<"Translation: "<<std::endl<<translation.transpose()<<reset<<std::endl;

  Eigen::Quaterniond q(rotation);
  std::cerr<<yellow<<"Quaternion: "<<q.coeffs().transpose()<<reset<<std::endl;

  // We use the points that is used for the CWC wrench cone. 
  for (const auto & point : getLocalContactAreaVerticies_())
  //for (const auto & point : pts)
  {
    // According to Arnaud, we have to use the inverse! 
    
    //inertialContactAreaVertices_.emplace_back( rotation* point + translation);

    inertialContactAreaVertices_.emplace_back( rotation * point + translation);
    //inertialContactAreaVertices_.emplace_back( rotation * point.translation() + translation);
    
    //inertialContactAreaVertices_.emplace_back( X_s_0.rotation() * point + X_s_0.translation());
    //inertialContactAreaVertices_.emplace_back( X_0_s.rotation().transpose() * point + X_0_s.translation());
  }


  // (1) Update the surface points: 
  inertialSurfaceVertices_.clear();
  //std::vector<Eigen::Vector3d> points;
  //const auto & surface = const_cast<mc_rbdyn::Surface &>(*contact.r1Surface());
  //const auto & surface = robot.surface(getContactParams().surfaceName);
    // Points in body frame
  //const auto & pts = surface.points();
  const sva::PTransformd & X_0_b = robot.bodyPosW(surface.bodyName());
  

  for(const auto & p : pts)
  {
    sva::PTransformd X_0_p = p * X_0_b;
    inertialSurfaceVertices_.push_back(X_0_p.translation());
    //std::cout<<green<<"point: "<<X_0_p.translation().transpose()<<reset<<std::endl;
  }

  // (2) Update the CoP 
  auto inputWrench= robot.forceSensor(getContactParams().sensorName).wrench();
  //std::cerr<<red<<"Void update of contact CoP "<<reset<<std::endl;
  
  measuredCoP_.setZero();
  if(inputWrench.force().z() > 10.0) 
  {
    // Don't calculate the CoP when the contact is not set. 
    measuredCoP_.x() = -inputWrench.couple().y() / inputWrench.force().z();
    measuredCoP_.y() = -inputWrench.couple().x() / inputWrench.force().z();
  }

   // Transform the CoP from the local frame to the surface frame.

  //measuredCoP_ = X_0_s.rotation().transpose()*measuredCoP_ + X_0_s.translation();
  measuredCoP_ = rotation*measuredCoP_ + translation;

 }

void McContact::update(const mc_rbdyn::Robot & robot)
{

  // (0) Update the CWC in the inertial frame 
  updateCWCInertial_(robot);

  // (1) Update the GraspMatrix 
  if(getContactParams().useSpatialVectorAlgebra)
  {
    calcGraspMatrix_(robot);
  }
  else
  {
    calcGeometricGraspMatrix_(robot);
  }

  // (2) Update the vertices 
  // In the future we may use GEOS to create more complex contact area. 
  updateContactAreaVerticiesAndCoP_(robot);


    
  // 
}

void McContact::updateCWCInertial_(const mc_rbdyn::Robot & robot)
{

  const sva::PTransformd & X_0_s = robot.surfacePose(getContactParams().surfaceName);

  auto rotationTranspose = X_0_s.rotation();
  auto translation = X_0_s.translation();

  // Update the Contact Wrench Cone (CWC) and the resultant wrench multiplier
  Eigen::Matrix6d rotationCorrection;
  rotationCorrection.setIdentity();


  rotationCorrection.block<3,3>(0,0) = rotationTranspose;
  rotationCorrection.block<3,3>(3,3) = rotationTranspose;

  CWCInertial_ = CWC_*rotationCorrection;

  resultantWrenchMultiplier_.setIdentity();
  resultantWrenchMultiplier_.block<3, 3>(3, 0) = crossMatrix(translation);


}

const McContact & McContactSet::getContact(const std::string & name) const
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

bool McContactSet::addContact(const McContactParams & inputParams, const mc_rbdyn::Robot & robot)
{
  auto opt = contacts_.find(inputParams.surfaceName);


  if(opt != (contacts_.end()))
  {
    std::cout << "Endeffector " << inputParams.surfaceName << " already exists." << std::endl;
    return false;
  }


  contacts_.insert({inputParams.surfaceName, {inputParams, robot}});

  //std::cout << "McContactSet: Adding end-effector contact: " << inputParams.surfaceName << std::endl;

  return true;
}

 

void McContactSet::update(const mc_rbdyn::Robot & robot)
{
  for(auto & contactPair:contacts_)
  {
    contactPair.second.update(robot);
  }
}

} // namespace mc_impact
