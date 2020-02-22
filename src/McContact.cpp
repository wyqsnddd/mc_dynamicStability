#include "McDynamicStability/McContact.h"

namespace mc_impact
{

McContact::McContact(const McContactParams & inputParams, const mc_rbdyn::Robot & robot) : mcContactParams_(inputParams), robot_(robot)
{
  // Initialize the local contact points with a square:

  localContactAreaVertices_.emplace_back(getContactParams().halfX, getContactParams().halfY, 0.0);
  localContactAreaVertices_.emplace_back(getContactParams().halfX, -getContactParams().halfY, 0.0);
  localContactAreaVertices_.emplace_back(-getContactParams().halfX, -getContactParams().halfY, 0.0);
  localContactAreaVertices_.emplace_back(-getContactParams().halfX, getContactParams().halfY, 0.0);

  initializeCWC_();

  resultantWrenchMultiplier_.setIdentity();

  if(getContactParams().initialContactStatus){
     setContact();
  } else {
     breakContact(); 
  }

#ifdef DEBUG
  std::cout << cyan << "Creating contact: " << getContactParams().surfaceName << reset << std::endl;
#endif

  update();
}

void McContact::calcGeometricGraspMatrix_()
{

  graspMatrix_.setIdentity();
  auto X_0_c = robot().surfacePose(getContactParams().surfaceName);

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

void McContact::calcGraspMatrix_()
{

  // --------------- Use the dual-matrix as the
  // sva::PTransformd X_c_0 = robot().surfacePose(getContactParams().surfaceName).inv();
  // graspMatrix_ = X_c_0.dualMatrix();

  // ---------------
  sva::PTransformd X_0_c = robot().surfacePose(getContactParams().surfaceName);
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

void McContact::updateContactAreaVerticiesAndCoP_()
{

  inertialContactAreaVertices_.clear();

  
  const sva::PTransformd & X_0_s = robot().surfacePose(getContactParams().surfaceName);
  //const sva::PTransformd & X_s_0 = robot().surfacePose(getContactParams().surfaceName).inv();

  // (0) Update the contact area points: 
  const auto & surface = robot().surface(getContactParams().surfaceName);
  const auto & pts = surface.points();
  auto rotation = X_0_s.rotation().inverse();
  //auto translation = - rotation * X_0_s.translation();
  auto translation = X_0_s.translation();




  /*
  std::cerr<<red<<"Updating contact: "<<getContactParams().surfaceName<<reset<<std::endl;
  std::cerr<<cyan<<"The surface: "<<surface.name()<<" has bodyName: "<<surface.bodyName()<<reset<<std::endl;

  std::cerr<<yellow<<"Rotation: "<<std::endl<<rotation<<reset<<std::endl;
  std::cerr<<yellow<<"Translation: "<<std::endl<<translation.transpose()<<reset<<std::endl;

  Eigen::Quaterniond q(rotation);
  std::cerr<<yellow<<"Quaternion: "<<q.coeffs().transpose()<<reset<<std::endl;
  */

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
  //const auto & surface = robot().surface(getContactParams().surfaceName);
    // Points in body frame
  //const auto & pts = surface.points();
  const sva::PTransformd & X_0_b = robot().bodyPosW(surface.bodyName());
  

  for(const auto & p : pts)
  {
    sva::PTransformd X_0_p = p * X_0_b;
    inertialSurfaceVertices_.push_back(X_0_p.translation());
    //std::cout<<green<<"point: "<<X_0_p.translation().transpose()<<reset<<std::endl;
  }

  // (2) Update the CoP 
  auto inputWrench= robot().forceSensor(getContactParams().sensorName).wrench();
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

void McContact::update()
{

  // (0) Update the CWC in the inertial frame 
  updateCWCInertial_();

  // (1) Update the GraspMatrix 
  if(getContactParams().useSpatialVectorAlgebra)
  {
    calcGraspMatrix_();
  }
  else
  {
    calcGeometricGraspMatrix_();
  }

  // (2) Update the vertices 
  // In the future we may use GEOS to create more complex contact area. 
  updateContactAreaVerticiesAndCoP_();


    
  // 
}

void McContact::updateCWCInertial_()
{

  const sva::PTransformd & X_0_s = robot().surfacePose(getContactParams().surfaceName);

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


void McContact::addGuiItems(mc_control::fsm::Controller &ctl) const
{

  mc_rtc::gui::ArrowConfig surfaceNormalConfig({0., 0., 1.});
  surfaceNormalConfig.start_point_scale = 0.0;
  surfaceNormalConfig.end_point_scale = 0.0;

  mc_rtc::gui::ArrowConfig surfaceXConfig({1., 0., 0.});
  surfaceXConfig.start_point_scale = 0.0;
  surfaceXConfig.end_point_scale = 0.0;

  mc_rtc::gui::ArrowConfig surfaceYConfig({0., 1., 0.});
  surfaceYConfig.start_point_scale = 0.0;
  surfaceYConfig.end_point_scale = 0.0;

  mc_rtc::gui::ArrowConfig surfaceZConfig({0., 0., 1.});
  surfaceZConfig.start_point_scale = 0.0;
  surfaceZConfig.end_point_scale = 0.0;

  double arrowLengthScale = 0.2;

	// TODO: new name for the contact points
  ctl.gui()->addElement(
        {"McContacts"},
        mc_rtc::gui::Polygon(getContactParams().surfaceName + "ContactVertices", mc_rtc::gui::Color(0., 0., 1.),
                             [this]() -> const std::vector<Eigen::Vector3d> & {
                                   return getContactAreaVertices();
                             }),
        mc_rtc::gui::Polygon(getContactParams().surfaceName + "SurfaceVertices", mc_rtc::gui::Color(0., 1.0, 0.0),
                             [this ]() -> const std::vector<Eigen::Vector3d> & {
                               return getSurfaceVertices();
                             }),
        mc_rtc::gui::Point3D(getContactParams().surfaceName + "CoP",
                             [this]() -> Eigen::Vector3d {
                               return measuredCop();
                             }),

        mc_rtc::gui::Arrow(getContactParams().surfaceName + "SurfaceNomal", surfaceNormalConfig,
                           [this]() -> Eigen::Vector3d {
                             // start of the arrow
                             auto X_0_s = this->robot().surfacePose(getContactParams().surfaceName);
                             return X_0_s.translation();
                           },
                           [this, arrowLengthScale]() {
                             // End of the arrow
                             auto X_0_s = this->robot().surfacePose(getContactParams().surfaceName);
                             Eigen::Vector3d z_axis =
                                 X_0_s.translation() + arrowLengthScale * (X_0_s.rotation().transpose().block<3, 1>(0, 2));
                             return z_axis;
                           }),
        mc_rtc::gui::Arrow(getContactParams().surfaceName + "SurfaceX", surfaceXConfig,
                           [this]() -> Eigen::Vector3d {
                             // start of the arrow
                             auto X_0_s = this->robot().surfacePose(getContactParams().surfaceName);
                             return X_0_s.translation();
                           },
                           [this, arrowLengthScale]() -> Eigen::Vector3d {
                             // End of the arrow
                             auto X_0_s = this->robot().surfacePose(getContactParams().surfaceName);
                             Eigen::Vector3d x_axis =
                                 X_0_s.translation() + arrowLengthScale * (X_0_s.rotation().transpose().block<3, 1>(0, 0));
                             return x_axis;
                           }),
        mc_rtc::gui::Arrow(getContactParams().surfaceName + "SurfaceY", surfaceYConfig,
                           [this]() -> Eigen::Vector3d {
                             // start of the arrow
                             auto X_0_s = this->robot().surfacePose(getContactParams().surfaceName);
                             return X_0_s.translation();
                           },
                           [this, arrowLengthScale]() -> Eigen::Vector3d {
                             // End of the arrow
                             auto X_0_s = this->robot().surfacePose(getContactParams().surfaceName);
                             Eigen::Vector3d y_axis =
                                 X_0_s.translation() + arrowLengthScale * (X_0_s.rotation().transpose().block<3, 1>(0, 1));
                             return y_axis;
                           }),
       mc_rtc::gui::Arrow(getContactParams().surfaceName + "SurfaceZ", surfaceZConfig,
                           [this]() -> Eigen::Vector3d {
                             // start of the arrow
                             auto X_0_s = this->robot().surfacePose(getContactParams().surfaceName);
                             return X_0_s.translation();
                           },
                           [this, arrowLengthScale]() -> Eigen::Vector3d {
                             // End of the arrow
                             auto X_0_s = this->robot().surfacePose(getContactParams().surfaceName);
                             Eigen::Vector3d y_axis =
                                 X_0_s.translation() + arrowLengthScale * (X_0_s.rotation().transpose().block<3, 1>(0, 2));
                             return y_axis;
                           })


    );
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

 

void McContactSet::update()
{
  for(auto & contactPair:contacts_)
  {
    contactPair.second.update();
  }
}

void McContactSet::addGuiItems(mc_control::fsm::Controller &ctl) const
{
  for(auto & contactPair:contacts_)
  {
    contactPair.second.addGuiItems(ctl);
  }
}

} // namespace mc_impact
