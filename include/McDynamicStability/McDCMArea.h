#pragma once

#include <mc_rbdyn/Robots.h>

#include "McDynamicStability/McContact.h"
#include "McDynamicStability/McPolytopeDescriptor.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <eigen-lssol/LSSOL_QP.h>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <polytope/staticStabilityPolytope.h>

#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/CoordinateSequenceFactory.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/Polygon.h>


#include "McComArea.h" 
#include "McZMPArea.h" 

namespace mc_impact
{

class McDCMArea
/*! \brief We calculate the `Multi-contact DCM-area` using the intersection of the `Multi-contat ZMP are` and the `Multi-contact Com static equilibrium` area
*/
{
  public:   
   McDCMArea(std::shared_ptr<mc_impact::McZMPArea<Eigen::Vector2d>> mcZMPAreaPtr,
                 std::shared_ptr<mc_impact::McComArea> mcComAreaPtr
);
   ~McDCMArea(){}

   /*! \brief It uses the internal pointers to the McComArea and the McZMPArea to calculate the McDCM area. 
    * We can use getPolygonVertices() to access the updated vertices belong to the intersection. 
    * \return false if there is no intersection between the ZMP area and the Com area.
    */
   bool updateMcDCMArea();

   /*! \return vertices of the McDCMArea.
    */
   inline const std::vector<Eigen::Vector2d> & getPolygonVertices() const
  {
    return polygonVertices_;
  }

  private:
  std::shared_ptr<McZMPArea<Eigen::Vector2d>> mcZMPAreaPtr_;
  std::shared_ptr<McComArea> mcComAreaPtr_;

  std::vector<Eigen::Vector2d> polygonVertices_;
};

}
