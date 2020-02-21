#pragma once

#include <polytope/staticStabilityPolytope.h>
#include <problemDescriptor/problemDescriptor.h>

namespace mc_impact
{

class McPolytopeDescriptor : public ProblemDescriptor
/*! \brief Inherite the interface class from "problemDescriptor" from stabiliPlus library
 */
{
public:
  McPolytopeDescriptor(std::string name = "Multi-contactPolytopeDescriptor") : ProblemDescriptor(name) {}
  ~McPolytopeDescriptor() {}

  // ----------- main class methods ----------
  /*! \brief Update the building blocks
   *
   */

  void update() override {}

  inline Eigen::MatrixXd & getA()
  {
    return m_A;
  }

  inline Eigen::VectorXd & getB()
  {
    return m_b;
  }
  inline Eigen::MatrixXd & getF()
  {
    return m_F;
  }
  inline Eigen::VectorXd & getf()
  {
    return m_f;
  }
};
} // namespace mc_impact
