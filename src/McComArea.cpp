#include "McDynamicStability/McComArea.h"

#include "McDynamicStability/Utils.h"

namespace mc_impact
{
// Repeat static constexpr declarations
// See also https://stackoverflow.com/q/8016780
constexpr double McComArea::LOWER_SLOPE;

constexpr double McComArea::UPPER_SLOPE;

constexpr int McComArea::GM_SIZE;

constexpr int McComArea::RM_SIZE;

McComArea::McComArea(const mc_rbdyn::Robot & robot,
                     std::shared_ptr<McContactSet> contactSetPtr,
                     const McComAreaParams & McComAreaParams)
: robot_(robot), contactsPtr_(contactSetPtr), McComAreaParams_(McComAreaParams)
{

  pdPtr_.reset(new McPolytopeDescriptor());

  // Fill in the initial values:
  polygonVertices_.emplace_back(-0.15, 0.12);
  polygonVertices_.emplace_back(-0.15, -0.12);
  polygonVertices_.emplace_back(0.15, -0.12);
  polygonVertices_.emplace_back(0.15, 0.12);

  // Initialize the vairables.
  updateMcComArea();

  std::cout << red << "McComArea is created." << reset << std::endl;
}

void McComArea::updateMcComArea()
{

  // Update the contacts (the grasp matrices will be updated):
  contactsPtr_->update(getRobot());

  // Update the Multi-contact Com static equilibrium area:
  computeMcComArea_();
}

void McComArea::computeMcComArea_()
{
  assert(getContactSet()->getContactMap().size() > 0);

  int numContact = static_cast<int>(getContactSet()->getContactMap().size());
  int rowCWC = static_cast<int>(getContactSet()->getContactMap().begin()->second.contactWrenchCone().rows());
  int colCWC = static_cast<int>(getContactSet()->getContactMap().begin()->second.contactWrenchCone().cols());

  ///-------------Part One: Equality constraint: C X = d
  Eigen::MatrixXd G;
  G.resize(GM_SIZE, GM_SIZE * numContact);
  G.setZero();

  // 1.1 We use this pair of matrix and vector to specify the MCWC constraint
  pdPtr_->getF().resize(numContact * rowCWC + 4, numContact * colCWC + 2);
  pdPtr_->getF().setZero();

  pdPtr_->getf().resize(numContact * rowCWC + 4);
  pdPtr_->getf().setZero();

  int count = 0;
  for(auto & contactPair : getContactSet()->getContactMap())
  {
    pdPtr_->getF().block(count * rowCWC, count * colCWC, rowCWC, colCWC) =
        contactPair.second.contactWrenchCone() * contactPair.second.getGraspMatrix();
    G.block<GM_SIZE, GM_SIZE>(0, count * GM_SIZE).setIdentity();

    count++;
  }

  // 1.2 These lines are added to limit the maximum search range.
  pdPtr_->getF().bottomRightCorner(4, 2)(0, 0) = 1.0;
  pdPtr_->getF().bottomRightCorner(4, 2)(1, 0) = -1.0;
  pdPtr_->getF().bottomRightCorner(4, 2)(2, 1) = 1.0;
  pdPtr_->getF().bottomRightCorner(4, 2)(3, 1) = -1.0;

  pdPtr_->getf().tail(4).setOnes();
  pdPtr_->getf().tail(4) = pdPtr_->getf().tail(4) * 100000;

  ///-------------Part Two: Equality constraint: C X = d
  int assumptionSize = 4;
  // Update the matrices according to the LIPM assumptions:
  /* (1) Sum of the external forces = -mg
   * (2) n dot torque = 0
   */
  if(!getParams().useLIPMAssumptions)
  {
    assumptionSize = 0;
  }

  pdPtr_->getA().resize(assumptionSize + 2, GM_SIZE * numContact + 2);
  pdPtr_->getA().setZero();
  pdPtr_->getB().resize(assumptionSize + 2);
  pdPtr_->getB().setZero();

  if(getParams().useLIPMAssumptions)
  {
    updateLIPMAssumptions_(numContact, G);
  }

  ///-------------Part Three: Projection E X = f
  /*
   * The equality constraints that specify the projection: polytope -> polygon (on the surface).
   * n \times torque = 0
   */

  // This is the Matrix E:

  Eigen::Matrix3d crossUz = crossMatrix(Eigen::Vector3d::UnitZ());
  // clang-format off
  pdPtr_->getA().block(assumptionSize, 0, 2, GM_SIZE * numContact) =
      1.0 / (getRobot().mass()* 9.81) 
      * crossUz.block<2, 3>(0, 0)
      * G.block(3, 0, 3, GM_SIZE * numContact);

  // clang-format on 
  
  // Move the varialbes to the same side:  
  pdPtr_->getA().block<2, 2>(assumptionSize, GM_SIZE * numContact) = -Eigen::Matrix2d::Identity();

  // This is the Vector f:
  //pdPtr_->getB().tail<2>() << 0.0, 0.0;


  ///-------------Part Four: projection 
  polytopeProjectorPtr_ = std::make_shared<StaticStabilityPolytope>(pdPtr_, getParams().iterationLimit, getParams().convergeThreshold, GLPK);

  polytopeProjectorPtr_->initSolver();

  polytopeProjectorPtr_->projectionStabilityPolyhedron();


  polygonVertices_.clear();

  if(getParams().debug)
  {
    std::cerr<< red<<"The stabiliPlus vertices are: " <<reset<< std::endl;
  
    for(auto & point : polytopeProjectorPtr_->getInnerVertices())
    {
	    std::cerr<<green<< point.transpose() << reset<<std::endl;
    }
    std::cerr<< red<< "End of StabiliPlus verticies." << reset<<std::endl;

  }
}


void McComArea::updateLIPMAssumptions_(int numContact, const Eigen::MatrixXd & inputG)
{

  // Set the vector d:
  pdPtr_->getB()(3) = - getRobot().mass()*9.81;

  Eigen::MatrixXd B;
  B.resize(4, 6);
  B.setZero();

  B.block<3, 3>(0, 0).setIdentity();
  B.block<1, 3>(3, 3) = Eigen::Vector3d::UnitZ().transpose();

  pdPtr_->getA().block(0, 0, 4, GM_SIZE * numContact) = B * inputG;

  if(getParams().debug)
  {
    std::cerr<< "Intermediate: Matrix  B is: " << std::endl << B << std::endl;
  }

}
};
