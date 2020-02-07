#include "McDynamicStability/McDCMArea.h"

namespace mc_impact
{
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


  // (1.2) Com Area:
  std::vector<Eigen::Vector2d> ComAreaVertices = mcComAreaPtr_->getPolygonVertices();  
  auto comPoints2dseq = factory.getCoordinateSequenceFactory()->create(static_cast<size_t>(0), 0);
  std::vector<geos::geom::Coordinate> comPoints;
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
  auto dcmSurfGeom = zmpSurfPoly->intersection(comSurfPoly.get());
  auto dcmSurfPoly = dynamic_cast<geos::geom::Polygon *>(dcmSurfGeom.get());
#else
  auto comSurfPoly = factory.createPolygon(std::move(comPoints2dshell), nullptr);
  auto dcmSurfGeom = zmpSurfPoly->intersection(comSurfPoly);
  auto dcmSurfPoly = dynamic_cast<geos::geom::Polygon *>(dcmSurfGeom);
#endif

  if(dcmSurfPoly == 0)
  {
    LOG_INFO( "Multi-contact ZMP area and COM static equilibrium area don't intersect.");
    return false; 
  }

  // (3) Output the vertices
  //std::vector<sva::PTransformd> res;
  auto dcmPoints = dcmSurfPoly->getExteriorRing()->getCoordinates();
  for(size_t i = 0; i < dcmPoints->getSize() - 1; ++i)
  {
    const geos::geom::Coordinate & p = dcmPoints->getAt(i);
    polygonVertices_.push_back({p.x, p.y});
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
  return true; 
}

}
