#include "McDynamicStability/McDCMArea.h"

namespace mc_impact
{
 constexpr double McDCMArea::LOWER_SLOPE;

 constexpr double McDCMArea::UPPER_SLOPE;



  McDCMArea::McDCMArea(std::shared_ptr<mc_impact::McZMPArea<Eigen::Vector2d>> mcZMPAreaPtr,
                 std::shared_ptr<mc_impact::McComArea> mcComAreaPtr):mcZMPAreaPtr_(mcZMPAreaPtr), mcComAreaPtr_(mcComAreaPtr)
  {
   // Fill in the initial values:
  polygonVertices_.emplace_back(-0.15, 0.12);
  polygonVertices_.emplace_back(-0.15, -0.12);
  polygonVertices_.emplace_back(0.15, -0.12);
  polygonVertices_.emplace_back(0.15, 0.12);

  updateMcDCMArea();
 
  std::cout << red << "McDCMArea is created." << reset << std::endl;
  }



bool McDCMArea::updateMcDCMArea()
{
  std::vector<Eigen::Vector2d> ZMPAreaVertices = mcZMPAreaPtr_->getPolygonVertices();  


  // (1) Create polygons with GEOS 
  
  // (1.1) ZMP polygon
  const geos::geom::GeometryFactory * factory_ptr = geos::geom::GeometryFactory::getDefaultInstance();
  const geos::geom::GeometryFactory & factory = *factory_ptr;

  auto zmpPoints2dseq = factory.getCoordinateSequenceFactory()->create(static_cast<size_t>(0), 0);
  std::vector<geos::geom::Coordinate> zmpPoints;
  zmpPoints.clear();
  for(const auto  & p : ZMPAreaVertices)
  {
     zmpPoints.push_back(geos::geom::Coordinate(p.x(), p.y()));
  }
  zmpPoints.push_back(geos::geom::Coordinate(zmpPoints[0]));
  zmpPoints2dseq->setPoints(zmpPoints);
  auto zmpPoints2dshell = factory.createLinearRing(std::move(zmpPoints2dseq));
#if GEOS_VERSION_MAJOR >= 3 && GEOS_VERSION_MINOR >= 8
  auto zmpSurfPoly = factory.createPolygon(std::move(zmpPoints2dshell));
#else
  auto zmpSurfPoly = factory.createPolygon(std::move(zmpPoints2dshell), nullptr);
#endif
  //std::cout<< red << "GEOS: Constructed ZMP polygon, which is valid? "<<zmpSurfPoly->isValid() << std::endl;
  
  auto zmpSurfPoly_two = zmpSurfPoly->buffer(0);
  /*
  if(zmpSurfPoly->isValid()) 
  {
    return false ;
  }
  */
  // (1.2) Com Area:
  std::vector<Eigen::Vector2d> ComAreaVertices = mcComAreaPtr_->getPolygonVertices();  
  auto comPoints2dseq = factory.getCoordinateSequenceFactory()->create(static_cast<size_t>(0), 0);
  std::vector<geos::geom::Coordinate> comPoints;
  comPoints.clear();
  for(const auto  & p : ComAreaVertices)
  {
    comPoints.push_back(geos::geom::Coordinate(p.x(), p.y()));
  }
  comPoints.push_back(geos::geom::Coordinate(comPoints[0]));
  comPoints2dseq->setPoints(comPoints);
  auto comPoints2dshell = factory.createLinearRing(std::move(comPoints2dseq));



  // (2) Calculate the intersection
#if GEOS_VERSION_MAJOR >= 3 && GEOS_VERSION_MINOR >= 8
  auto comSurfPoly = factory.createPolygon(std::move(comPoints2dshell));
#else
  auto comSurfPoly = factory.createPolygon(std::move(comPoints2dshell), nullptr);
#endif
  auto comSurfPoly_two = comSurfPoly->buffer(0);
  /*
  if(comSurfPoly->isValid()) 
  {
    return false ;
  }
  */

 //std::cout<< red << "GEOS: Constructed COM polygon, which is valid? "<< comSurfPoly->isValid() <<reset <<std::endl;

#if GEOS_VERSION_MAJOR >= 3 && GEOS_VERSION_MINOR >= 8
  auto dcmSurfGeom = comSurfPoly_two->intersection(zmpSurfPoly_two.get());
  auto dcmSurfPoly = dynamic_cast<geos::geom::Polygon *>(dcmSurfGeom.get());
#else
  auto dcmSurfGeom = comSurfPoly_two->intersection(zmpSurfPoly_two);
  auto dcmSurfPoly = dynamic_cast<geos::geom::Polygon *>(dcmSurfGeom);
#endif



  if(mcZMPAreaPtr_->getParams().debug)
  {
    std::cout<< red << "GEOS: Computed intersection "<<reset <<std::endl;
  }
  if(dcmSurfPoly == 0)
  {
    LOG_INFO( "Multi-contact ZMP area and COM static equilibrium area don't intersect.");
    return false; 
  }

  // (3) Output the vertices
  //std::vector<sva::PTransformd> res;
  auto dcmPoints = dcmSurfPoly->getExteriorRing()->getCoordinates();
  polygonVertices_.clear();
  for(size_t i = 0; i < dcmPoints->getSize() - 1; ++i)
  {
    const geos::geom::Coordinate & p = dcmPoints->getAt(i);
    polygonVertices_.push_back({p.x, p.y});
  }
  if(mcZMPAreaPtr_->getParams().debug)
  {
    std::cout<< red << "GEOS: extracted DCM vertices. "<<reset <<std::endl;
  }

  if(mcZMPAreaPtr_->getParams().debug)
  {
    std::cerr << cyan << "The DCM-area vertices are: " << reset << std::endl;
    for(auto & point : getPolygonVertices())
    {
      std::cerr << point.transpose() << std::endl;
    }
    std::cerr << "--------------------------------" << std::endl;
  }
 
  pointsToInequalityMatrix<Eigen::Vector2d>(polygonVertices_, ieqConstraintBlocks_.G, ieqConstraintBlocks_.h, LOWER_SLOPE, UPPER_SLOPE);

  return true; 
}

}
