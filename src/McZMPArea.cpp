#include "McDynamicStability/McZMPArea.h"

#include "McDynamicStability/Utils.h"

namespace mc_impact
{

// Repeat static constexpr declarations
// See also https://stackoverflow.com/q/8016780
template<typename Point>
constexpr double McZMPArea<Point>::LOWER_SLOPE;

template<typename Point>
constexpr double McZMPArea<Point>::UPPER_SLOPE;

template<typename Point>
constexpr int McZMPArea<Point>::GM_SIZE;

template<typename Point>
constexpr int McZMPArea<Point>::RM_SIZE;

template<typename Point>
McZMPArea<Point>::McZMPArea(const mc_rbdyn::Robot & robot,
                            std::shared_ptr<McContactSet> contactSetPtr,
                            const McZMPAreaParams & mcZMPAreaParams)
: robot_(robot), contactsPtr_(contactSetPtr), McZMPAreaParams_(mcZMPAreaParams)
{
  // pdPtr_ = std::make_shared<McPolytopeDescriptor>();
  pdPtr_.reset(new McPolytopeDescriptor());

  /*
  polytopeProjectorPtr_ = std::make_shared<StaticStabilityPolytope>(pdPtr_, getParams().iterationLimit,
                                                                    getParams().convergeThreshold, GLPK);

								    */
  // Fill in the initial values:
  /*
  polygonVertices_.emplace_back(-0.10, 0.18);
  polygonVertices_.emplace_back(0.14, 0.18);
  polygonVertices_.emplace_back(0.14, -0.18);
  polygonVertices_.emplace_back(-0.10, -0.18);
  */

  /*
  polygonVertices_.emplace_back(-0.12, 0.1);
  polygonVertices_.emplace_back(-0.12, -0.1);
  polygonVertices_.emplace_back(0.12, -0.1);
  polygonVertices_.emplace_back(0.12, 0.1);
  */

  polygonVertices_.emplace_back(-0.15, 0.12);
  polygonVertices_.emplace_back(-0.15, -0.12);
  polygonVertices_.emplace_back(0.15, -0.12);
  polygonVertices_.emplace_back(0.15, 0.12);



  // Initialize the vairables.
  computeMcZMPArea();

  std::cout << "McZMPArea is created." << std::endl;
}

template<typename Point>
void McZMPArea<Point>::computeMcZMPArea(double height)
{

  assert(getContactSet()->getContactMap().size() > 0);

  int numContact = static_cast<int>(getContactSet()->getContactMap().size());
  int rowCWC = static_cast<int>(getContactSet()->getContactMap().begin()->second.contactWrenchCone().rows());
  int colCWC = static_cast<int>(getContactSet()->getContactMap().begin()->second.contactWrenchCone().cols());

  double mass = getRobot().mass(); // The result does not depend on the mass.
  Eigen::MatrixXd G;
  G.resize(GM_SIZE, GM_SIZE * numContact);
  G.setZero();

  // Polytope:
  // Eigen::MatrixXd F; ///< Inequality constraint of the Halfspace representation of the polytope
  // Eigen::VectorXd b;

  // Eigen::MatrixXd & F = pdPtr_->getA();

  pdPtr_->getF().resize(numContact * rowCWC + 4, numContact * colCWC + 2);
  pdPtr_->getF().setZero();

  // Eigen::VectorXd & b = pdPtr_->getB();
  // b.resize(numContact * colCWC);
  pdPtr_->getf().resize(numContact * rowCWC + 4);
  pdPtr_->getf().setZero();

  int count = 0;
  for(auto & contactPair : getContactSet()->getContactMap())
  {
    Eigen::Matrix6d tempGraspMatrix;
    contactPair.second.calcGeometricGraspMatrix(tempGraspMatrix, getRobot());
    //contactPair.second.calcGraspMatrix(tempGraspMatrix, getRobot());
    //std::cout<<"The local cwc matrix is: "<< std::endl<<contactPair.second.contactWrenchCone()<<std::endl;
    //std::cout<<"The grasp matrix is: "<< std::endl<< tempGraspMatrix <<std::endl;

    pdPtr_->getF().block(count * rowCWC, count * colCWC, rowCWC, colCWC) =
        contactPair.second.contactWrenchCone() * tempGraspMatrix;
    //std::cout<<"The cwc matrix is: "<< std::endl<<pdPtr_->getF().block(count * rowCWC, count * colCWC, rowCWC, colCWC)<<std::endl;

    //G.block<GM_SIZE, GM_SIZE>(0, count * GM_SIZE) = tempGraspMatrix;
    G.block<GM_SIZE, GM_SIZE>(0, count * GM_SIZE).setIdentity();
    count++;
  }

  pdPtr_->getF().bottomRightCorner(4,2)(0,0) = 1.0;
  pdPtr_->getF().bottomRightCorner(4,2)(1,0) = -1.0;
  pdPtr_->getF().bottomRightCorner(4,2)(2,1) = 1.0;
  pdPtr_->getF().bottomRightCorner(4,2)(3,1) = -1.0;
	  
  pdPtr_->getf().tail(4).setOnes();
  pdPtr_->getf().tail(4) = pdPtr_->getf().tail(4)*100000;

  ///-------------Part: Equality constraint: C X = d

  // Eigen::MatrixXd C; ///< Equality constraint that specifies the assumptions: (1) z acceleration = 0.0 (2) zero
  // angular momentum.
  // Eigen::Vector4d d;
 

  int assumptionSize = 4;
    
  // Update the matrices according to the LIPM assumptions: 
  if(! getParams().useLIPMAssumptions)
  {
    assumptionSize = 0;	  

  }

  pdPtr_->getA().resize(assumptionSize + 2, GM_SIZE * numContact + 2);
  pdPtr_->getA().setZero();
  pdPtr_->getB().resize(assumptionSize + 2);
  pdPtr_->getB().setZero();

  if( getParams().useLIPMAssumptions)
  {
    updateLIPMAssumptions_(numContact, G);
  }



    ///-------------Part: Projection E X = f
  // Projection
  // Eigen::MatrixXd E; ///< The equality constraints that specify the projection: polytope -> polygon (on the surface).
  // Eigen::Vector2d f;

  // Matrix E:

  pdPtr_->getA().block(assumptionSize, 0, 2, GM_SIZE * numContact) =
      (height - getRobot().com().z()) / (mass * 9.81) * G.block(0, 0, 2, GM_SIZE * numContact);
  pdPtr_->getA().block<2, 2>(assumptionSize, GM_SIZE * numContact) = -Eigen::Matrix2d::Identity();

  /*
  pdPtr_->getA().block(0, 0, 2, GM_SIZE*numContact) = (height - getRobot().com().z()) / (mass * 9.81) * G.block(0, 0, 2,
  GM_SIZE*numContact); pdPtr_->getA().block<2,2>(0, GM_SIZE*numContact) = -Eigen::Matrix2d::Identity();
  */

  // Vector f:
  pdPtr_->getB().tail<2>() << -getRobot().com().x(), -getRobot().com().y();



  //------------------------ temporary fix
/*  
  pdPtr_->getF().resize(numContact * rowCWC + 4, numContact * colCWC + 2);
  pdPtr_->getF().setZero();
  pdPtr_->getf().resize(numContact * rowCWC + 4);
  pdPtr_->getf().setZero();



  pdPtr_->getF()<<
  1.00005630e+00, -8.35276168e-04, -4.94860272e-01,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  
 -9.99943513e-01, -1.65727323e-03, -4.95086077e-01,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
 -3.54317128e-04,  9.98750471e-01, -4.97491076e-01,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  4.67108777e-04, -1.00124302e+00, -4.92455273e-01,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  7.40589016e-06,  9.31449418e-04, -1.44648894e-01,  9.99999909e-01,  4.10998530e-04,  1.12902285e-04,  0.00000000e+00,  0.00000000e+00,
  7.40589016e-06, -1.25877058e-03,  1.46493065e-02, -9.99999909e-01, -4.10998530e-04, -1.12902285e-04,  0.00000000e+00,  0.00000000e+00,
 -1.08143759e-03, -3.02142614e-04, -1.03085019e-01, -4.10712952e-04,  9.99996746e-01, -2.51790171e-03,  0.00000000e+00,  0.00000000e+00,
  1.10878241e-03, -3.02142614e-04, -1.36914219e-01,  4.10712952e-04, -9.99996746e-01,  2.51790171e-03,  0.00000000e+00,  0.00000000e+00,
  1.45152294e-01,  1.02339112e-01, -6.08128536e-02, -4.94657473e-01, -4.97694425e-01, -9.98806410e-01,  0.00000000e+00,  0.00000000e+00,
  1.44166761e-01, -1.37660107e-01, -4.34639575e-02, -4.95064058e-01,  4.92251847e-01, -1.00129901e+00,  0.00000000e+00,  0.00000000e+00,
  1.51523054e-02,  1.03369786e-01, -1.39676117e-01,  4.95291931e-01, -4.97287557e-01, -9.98694642e-01,  0.00000000e+00,  0.00000000e+00,
  1.41667729e-02, -1.36629433e-01, -1.22327221e-01,  4.94885346e-01,  4.92658715e-01, -1.00118724e+00,  0.00000000e+00,  0.00000000e+00,
 -1.44047324e-01, -1.03884337e-01, -6.02232345e-02, -4.94885346e-01, -4.92658715e-01,  1.00118724e+00,  0.00000000e+00,  0.00000000e+00,
 -1.45229998e-01,  1.36114882e-01, -4.40829312e-02, -4.95291931e-01,  4.97287557e-01,  9.98694642e-01,  0.00000000e+00,  0.00000000e+00,
 -1.40473353e-02, -1.02746804e-01, -1.39057143e-01,  4.95064058e-01, -4.92251847e-01,  1.00129901e+00,  0.00000000e+00,  0.00000000e+00,
 -1.52300100e-02,  1.37252415e-01, -1.22916840e-01,  4.94657473e-01,  4.97694425e-01,  9.98806410e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00, -1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  1.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00, -1.00000000e+00;


 pdPtr_->getf()<<
      0., 
      0.,
      0.,
      0.,
      0.,
      0.,
      0.,
      0.,
      0.,
      0.,
      0.,
      0.,
      0.,
      0.,
      0.,
      0.,
 100000.,
 100000.,
 100000.,
 100000.;

  pdPtr_->getA().resize(assumptionSize + 2, GM_SIZE * numContact + 2);
  pdPtr_->getA().setZero();
  pdPtr_->getB().resize(assumptionSize + 2);
  pdPtr_->getB().setZero();


  pdPtr_->getA()<<
 -1.89616533e-03, -7.78228838e-07,  4.08390802e-05, -9.96827709e-07,  2.42705875e-03, -6.11111527e-06,  0.00000000e+00,  0.00000000e+00, 
  7.77688096e-07, -1.89615934e-03,  1.98081340e-04, -2.42706643e-03, -9.97520823e-07, -2.74021371e-07,  0.00000000e+00,  0.00000000e+00,
 -2.15740143e-07,  4.76757778e-06,  1.89350158e-03,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
 -1.93202545e-04,  6.42846204e-05, -5.84817269e-08, -2.76532139e-07,  6.11100217e-06,  2.42705894e-03,  0.00000000e+00,  0.00000000e+00,
 -2.96062543e-03, -1.21681281e-06, -3.34261407e-07,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00, -1.00000000e+00,  0.00000000e+00,
  1.21596732e-06, -2.96061606e-03,  7.45456451e-06,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00, -1.00000000e+00;

  pdPtr_->getB()<<
  9.57196e-03, 
 -4.97186e-05,
  7.80163e-01,
  0.00000e+00,
 -9.57196e-03,
  4.97186e-05;
*/
  // pdPtr_->getB() << -getRobot().com().x(), -getRobot().com().y();

  std::cout << "A matrix is: " << std::endl << pdPtr_->getA() << std::endl;
  std::cout << "A matrix size is: " << std::endl << pdPtr_->getA().rows() << ", " << pdPtr_->getA().cols() << std::endl;
  std::cout << "B vector is: " << std::endl << pdPtr_->getB().transpose() << std::endl;
  std::cout << "B vector size is: " << std::endl << pdPtr_->getB().size() << std::endl;

  std::cout << "F matrix is: " << std::endl << pdPtr_->getF() << std::endl;
  std::cout << "F matrix size is: " << std::endl << pdPtr_->getF().rows() << ", " << pdPtr_->getF().cols() << std::endl;
  std::cout << "f vector is: " << std::endl << pdPtr_->getf().transpose() << std::endl;
  std::cout << "f vector size is: " << std::endl << pdPtr_->getf().size() << std::endl;

  std::cout << "Intermediate: Matrix  G is: " << std::endl << G << std::endl;

  polytopeProjectorPtr_ = std::make_shared<StaticStabilityPolytope>(pdPtr_, getParams().iterationLimit,
                                                                    getParams().convergeThreshold, GLPK);


  polytopeProjectorPtr_->initSolver();

  polytopeProjectorPtr_->projectionStabilityPolyhedron();

  numVertex_ = static_cast<int>(polytopeProjectorPtr_->getPolygonVerticies().size());
  // numVertex_ = static_cast<int>(polytopeProjectorPtr_->constraintPlanes().size());

  polygonVertices_.clear();

  std::cout << "The stabiliPlus vertices are: " << std::endl;
  for(auto & point : polytopeProjectorPtr_->getPolygonVerticies())
  {
    std::cout << point->innerVertex().transpose() << " with search direction:"<< point->searchDir().transpose()<<std::endl;
  }

  std::cout<<"End of StabiliPlus verticies."<<std::endl;


  // Construct the matrices:
  pointsToInequalityMatrix<StaticPoint>(polytopeProjectorPtr_->getPolygonVerticies(), ieqConstraintBlocks_.G_zmp,
                                        ieqConstraintBlocks_.h_zmp, polygonVertices_, LOWER_SLOPE, UPPER_SLOPE);
  std::cout << "The projected vertices are: " << std::endl;
  for(auto & point : polygonVertices_)
  {
    std::cout << point.transpose() << std::endl;
  }

  std::cout << "--------------------------------" << std::endl;
}

template<typename Point>
void McZMPArea<Point>::updateLIPMAssumptions_(int numContact, const Eigen::MatrixXd & inputG)
{
  pdPtr_->getB().head<3>() = getRobot().com();

  Eigen::Matrix3d crossUz = crossMatrix(Eigen::Vector3d::UnitZ());

  Eigen::MatrixXd B;
  B.resize(4, 6);
  B.setZero();
  B.block<3, 3>(0, 0).setIdentity();
  B.block<3, 3>(0, 0) *= getRobot().com().z();
  B.block<3, 3>(0, 3) = crossUz;
  //B.block<1, 3>(3, 3) = getRobot().com().transpose();
  
  /*
  B.block<1, 3>(3, 0) = -(crossUz * getRobot().com()).transpose();
  B.block<1, 3>(3, 3) = Eigen::Vector3d::UnitZ().transpose();
  */
 
  B.block<1, 3>(3, 3) = getRobot().com().transpose();

  // Matrix C:
  pdPtr_->getA().block(0, 0, 4, GM_SIZE * numContact) = 1.0 / (getRobot().mass()* 9.81) * (B * inputG);


  std::cout << "Intermediate: Matrix  B is: " << std::endl << B << std::endl;

}

// Instantiate the McZMPArea
template class mc_impact::McZMPArea<Eigen::Vector3d>;
template class mc_impact::McZMPArea<Eigen::Vector2d>;

} // namespace mc_impact
