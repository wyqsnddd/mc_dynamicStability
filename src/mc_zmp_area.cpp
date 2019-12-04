# include "mc_zmp_area.h"


namespace mc_impact{

template<typename Point>
McZMPArea<Point>::McZMPArea(const mc_rbdyn::Robot & robot, const struct McZMPAreaParams params):robot_(robot), params_(params)
{

  std::cout<<"McZMPArea is created."<<std::endl;
}

template<typename Point>
void McZMPArea<Point>::computeMcZMPArea(std::vector<Point> & zmpVerticies , double height )
{

  assert(contacts.getContactMap().size()>0);

  int numContact = static_cast<int>(contacts.getContactMap().size());
  int rowCWC = static_cast<int>(contacts.getContactMap().begin()->second.contactWrenchCone().rows());
  int colCWC = static_cast<int>(contacts.getContactMap().begin()->second.contactWrenchCone().cols());

  double mass = getRobot().mass(); // The result does not depend on the mass.

  Eigen::Matrix6d G;
  G.resize(6, 6*numContact);
  G.setZero();

  // Polytope: 
  Eigen::MatrixXd F; ///< Inequality constraint of the Halfspace representation of the polytope
  Eigen::VectorXd b;

  F.resize(numContact*rowCWC, numContact*colCWC);
  F.setZero();

  b.resize(numContact*colCWC);
  b.setZero();

  int count = 0; 
  for(auto & contactPair:contacts.getContactMap())
  {
    Eigen::Matrix6d rotationToWorld; 
    rotationToWorld.setIdentity();

    Eigen::Vector3d contactTranslation =  getRobot().surface(contactPair.second.getContactParams().surfaceName).X_0_s(getRobot()).translation();

    rotationToWorld.block<3,3>(0, 0) = getRobot().surface(contactPair.second.getContactParams().surfaceName).X_0_s(getRobot()).rotation();
 

    rotationToWorld.block<3,3>(3, 3) =  rotationToWorld.block<3,3>(0, 0);

    F.block(count*rowCWC, count*colCWC, rowCWC, colCWC) = contactPair.second.contactWrenchCone()*rotationToWorld;


    Eigen::Matrix6d graspMatrix;
    graspMatrix.setIdentity();
    graspMatrix.block<3,3>(3,0) =  crossMatrix_(contactTranslation - Eigen::Vector3d::Zero());
    graspMatrix *= rotationToWorld;

    G.block<6,6>(0, numContact*6) = graspMatrix; 

    count++;
  }
  
  Eigen::MatrixXd C; ///< Equality constraint that specifies the assumptions: (1) z acceleration = 0.0 (2) zero angular momentum.
  Eigen::Vector4d d;
  d.setZero();
  d.segment<3>(0) = getRobot().com();

  C.resize(4,6);
  Eigen::Matrix3d crossUz = crossMatrix_(Eigen::Vector3d::UnitZ());
 
  Eigen::MatrixXd B;
  B.resize(4,6);
  B.setZero();
  B.block<3,3>(0,0) *= getRobot().com().z();
  B.block<3,3>(0,3) =  crossUz; 
  B.block<1,3>(3,3) =  getRobot().com().transpose(); 
  B.block<1,3>(3,0) =  -(crossUz*getRobot().com()).transpose(); 

  C = 1.0/(mass*9.81)*(B*G);
 
  // Projection
  Eigen::MatrixXd E; ///< The equality constraints that specify the projection: polytope -> polygon (on the surface).
  Eigen::Vector2d f;

  E = (height - getRobot().com().z())/(mass*9.81)*G.block<2,6>(0,0);
  f << getRobot().com().x(), getRobot().com().y();


}

template<typename Point>
Eigen::Matrix3d McZMPArea<Point>::crossMatrix_(const Eigen::Vector3d &  input)
{

  Eigen::Matrix3d skewSymmetricMatrix = Eigen::Matrix3d::Zero();

  skewSymmetricMatrix(0,1) = -input(2);
  skewSymmetricMatrix(1,0) = input(2);

  skewSymmetricMatrix(0,2) = input(1);
  skewSymmetricMatrix(2,0) = -input(1);
  
 
  skewSymmetricMatrix(1,2) = -input(0);
  skewSymmetricMatrix(2,1) = input(0);
  
  return skewSymmetricMatrix;
}
// The explicit instantiation
template struct mc_impact::McZMPArea<Eigen::Vector3d>;
template struct mc_impact::McZMPArea<Eigen::Vector2d>;

}
