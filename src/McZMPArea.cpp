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

  polytopeProjectorPtr_ = std::make_shared<StaticStabilityPolytope>(pdPtr_, getParams().iterationLimit,
                                                                    getParams().convergeThreshold, GLPK);

  // Fill in the initial values:
  polygonVertices_.emplace_back(-0.10, 0.18);
  polygonVertices_.emplace_back(0.14, 0.18);
  polygonVertices_.emplace_back(0.14, -0.18);
  polygonVertices_.emplace_back(-0.10, -0.18);

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

  pdPtr_->getF().resize(numContact * rowCWC, numContact * colCWC + 2);
  pdPtr_->getF().setZero();

  // Eigen::VectorXd & b = pdPtr_->getB();
  // b.resize(numContact * colCWC);
  pdPtr_->getf().resize(numContact * rowCWC);
  pdPtr_->getf().setZero();

  int count = 0;
  for(auto & contactPair : getContactSet()->getContactMap())
  {
    Eigen::Matrix6d tempGraspMatrix;
    contactPair.second.calcGeometricGraspMatrix(tempGraspMatrix, getRobot());

    pdPtr_->getF().block(count * rowCWC, count * colCWC, rowCWC, colCWC) =
        contactPair.second.contactWrenchCone() * tempGraspMatrix;

    G.block<GM_SIZE, GM_SIZE>(0, count * RM_SIZE) = tempGraspMatrix;
    count++;
  }

  /*
    for(auto & contactPair : getContactSet()->getContactMap())
    {
      Eigen::Matrix6d rotationToWorld;
      rotationToWorld.setIdentity();

      Eigen::Vector3d contactTranslation =
          getRobot().surface(contactPair.second.getContactParams().surfaceName).X_0_s(getRobot()).translation();

      rotationToWorld.block<RM_SIZE, RM_SIZE>(0, 0) =
          getRobot().surface(contactPair.second.getContactParams().surfaceName).X_0_s(getRobot()).rotation().transpose();

      rotationToWorld.block<RM_SIZE, RM_SIZE>(3, 3) = rotationToWorld.block<RM_SIZE, RM_SIZE>(0, 0);

      pdPtr_->getF().block(count * rowCWC, count * colCWC, rowCWC, colCWC) = contactPair.second.contactWrenchCone() *
  rotationToWorld;

      Eigen::Matrix6d graspMatrix;

      Eigen::Vector3d testTranslation =
          getRobot().surface(contactPair.second.getContactParams().surfaceName).X_0_s(getRobot()).inv().translation();

      graspMatrix.setIdentity();
      graspMatrix.block<RM_SIZE, RM_SIZE>(3, 0) = crossMatrix(testTranslation);
      //graspMatrix.block<RM_SIZE, RM_SIZE>(3, 0) = crossMatrix_(contactTranslation - Eigen::Vector3d::Zero());

      //graspMatrix *= rotationToWorld;

      G.block<GM_SIZE, GM_SIZE>(0, count* RM_SIZE) = graspMatrix;
  count++;
    }
  */
  ///-------------Part: Equality constraint: C X = d

  // Eigen::MatrixXd C; ///< Equality constraint that specifies the assumptions: (1) z acceleration = 0.0 (2) zero
  // angular momentum.
  // Eigen::Vector4d d;

  pdPtr_->getA().resize(6, GM_SIZE * numContact + 2);
  pdPtr_->getA().setZero();
  pdPtr_->getB().resize(6);
  pdPtr_->getB().setZero();

  /*
  pdPtr_->getA().resize(2, GM_SIZE*numContact + 2);
  pdPtr_->getA().setZero();
  pdPtr_->getB().resize(2);
  pdPtr_->getB().setZero();

  */
  // Vector d
  pdPtr_->getB().head<3>() = getRobot().com();

  Eigen::Matrix3d crossUz = crossMatrix(Eigen::Vector3d::UnitZ());

  Eigen::MatrixXd B;
  B.resize(4, 6);
  B.setZero();
  B.block<3, 3>(0, 0).setIdentity();
  B.block<3, 3>(0, 0) *= getRobot().com().z();
  B.block<3, 3>(0, 3) = crossUz;
  B.block<1, 3>(3, 3) = getRobot().com().transpose();
  // B.block<1, 3>(3, 0) = -(crossUz * getRobot().com()).transpose();

  // Matrix C:
  pdPtr_->getA().block(0, 0, 4, GM_SIZE * numContact) = 1.0 / (mass * 9.81) * (B * G);

  ///-------------Part: Projection E X = f
  // Projection
  // Eigen::MatrixXd E; ///< The equality constraints that specify the projection: polytope -> polygon (on the surface).
  // Eigen::Vector2d f;

  // Matrix E:

  pdPtr_->getA().block(4, 0, 2, GM_SIZE * numContact) =
      (height - getRobot().com().z()) / (mass * 9.81) * G.block(0, 0, 2, GM_SIZE * numContact);
  pdPtr_->getA().block<2, 2>(4, GM_SIZE * numContact) = -Eigen::Matrix2d::Identity();

  /*
  pdPtr_->getA().block(0, 0, 2, GM_SIZE*numContact) = (height - getRobot().com().z()) / (mass * 9.81) * G.block(0, 0, 2,
  GM_SIZE*numContact); pdPtr_->getA().block<2,2>(0, GM_SIZE*numContact) = -Eigen::Matrix2d::Identity();
  */

  // Vector f:
  pdPtr_->getB().tail<2>() << -getRobot().com().x(), -getRobot().com().y();
  // pdPtr_->getB() << -getRobot().com().x(), -getRobot().com().y();

  std::cout << "A matrix is: " << std::endl << pdPtr_->getA() << std::endl;
  std::cout << "A matrix size is: " << std::endl << pdPtr_->getA().rows() << ", " << pdPtr_->getA().cols() << std::endl;
  std::cout << "B vector is: " << std::endl << pdPtr_->getB().transpose() << std::endl;
  std::cout << "B vector size is: " << std::endl << pdPtr_->getB().size() << std::endl;

  std::cout << "F matrix is: " << std::endl << pdPtr_->getF() << std::endl;
  std::cout << "F matrix size is: " << std::endl << pdPtr_->getF().rows() << ", " << pdPtr_->getF().cols() << std::endl;
  std::cout << "f vector is: " << std::endl << pdPtr_->getf().transpose() << std::endl;
  std::cout << "f vector size is: " << std::endl << pdPtr_->getf().size() << std::endl;

  std::cout << "Intermediate: Matrix  B is: " << std::endl << B << std::endl;
  std::cout << "Intermediate: Matrix  G is: " << std::endl << G << std::endl;

  polytopeProjectorPtr_->initSolver();

  polytopeProjectorPtr_->projectionStabilityPolyhedron();

  numVertex_ = static_cast<int>(polytopeProjectorPtr_->getPolygonVerticies().size());
  // numVertex_ = static_cast<int>(polytopeProjectorPtr_->constraintPlanes().size());

  polygonVertices_.clear();

  std::cout << "The stabiliPlus vertices are: " << std::endl;
  for(auto & point : polytopeProjectorPtr_->getPolygonVerticies())
  {
    std::cout << point->innerVertex().transpose() << std::endl;
  }

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

// Instantiate the McZMPArea
template class mc_impact::McZMPArea<Eigen::Vector3d>;
template class mc_impact::McZMPArea<Eigen::Vector2d>;

} // namespace mc_impact
