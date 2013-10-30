//------------------------------------------------------------------------
// BoostTessellator
// 
// Polytope wrapper for the native 2D Voronoi tessellator in Boost.Polygon
// v1.52 or greater
//------------------------------------------------------------------------
#ifndef __Polytope_BoostTessellator__
#define __Polytope_BoostTessellator__

#if HAVE_BOOST_VORONOI

#include <vector>
#include <cmath>

#include "Tessellator.hh"
#include "Clipper2d.hh"
#include "BoostOrphanage.hh"
#include "Point.hh"
#include "polytope_tessellator_utilities.hh"

// The Voronoi tools in Boost.Polygon
#include <boost/polygon/voronoi.hpp>

// We use the Boost.Geometry library to handle polygon intersections and such.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

// Some useful typedefs
typedef int64_t CoordHash;
typedef std::pair<int, int> EdgeHash;
typedef polytope::Point2<CoordHash> IntPoint;
typedef polytope::Point2<double> RealPoint;
typedef RealPoint PointType;

// Boost.Geometry typedefs
typedef boost::geometry::model::polygon<PointType, false> BGpolygon;
typedef boost::geometry::model::ring   <PointType, false> BGring;
typedef boost::geometry::model::polygon<IntPoint , false> IntPolygon;
typedef boost::geometry::model::ring   <IntPoint , false> IntRing;

//------------------------------------------------------------------------
// Map Polytope's point class to Boost.Polygon
//------------------------------------------------------------------------
namespace boost{
namespace polygon{

template <>
struct geometry_concept<IntPoint> { typedef point_concept type; };
  
template <>
struct point_traits<IntPoint> {
  typedef CoordHash coordinate_type;
   
  static inline coordinate_type get(const IntPoint& point, orientation_2d orient) {
    return (orient == HORIZONTAL) ? point.x : point.y;
  }
};


// //------------------------------------------------------------------------
// // Custom comparison operator
// //------------------------------------------------------------------------
// struct polytope_ulp_comparison {
//   enum Result {
//     LESS  = -1,
//     EQUAL =  0,
//     MORE  =  1
//   };
//    
//   Result operator()(CoordHash a, CoordHash b, unsigned int maxUlps) const {
//     if (a == b)
//       return EQUAL;
//     if (a >  b) {
//       Result res = operator()(b, a, maxUlps);
//       if (res == EQUAL) return res;
//       return (res == LESS) ? MORE : LESS;
//     }
//   }
// };
// 
// 
// //------------------------------------------------------------------------
// // Custom voronoi diagram traits
// //------------------------------------------------------------------------
// struct polytope_voronoi_diagram_traits {
//    typedef CoordHash coordinate_type;
//    typedef voronoi_cell<coordinate_type> cell_type;
//    typedef voronoi_vertex<coordinate_type> vertex_type;
//    typedef voronoi_edge<coordinate_type> edge_type;
//    typedef struct {
//    public:
//       enum {ULPS = 128};
//       bool operator()(const vertex_type &v1, const vertex_type &v2) const {
//          return (ulp_cmp(v1.x(), v2.x(), ULPS) == polytope_ulp_comparison::EQUAL and
//                  ulp_cmp(v1.y(), v2.y(), ULPS) == polytope_ulp_comparison::EQUAL);
//       }
//    private:
//       polytope_ulp_comparison ulp_cmp;
//    } vertex_equality_predicate_type;
// };

} //end boost namespace
} //end polygon namespace


//------------------------------------------------------------------------
// The Boost.Polygon Voronoi diagram object
//------------------------------------------------------------------------
// typedef boost::polygon::voronoi_diagram
//   <CoordHash,boost::polygon::polytope_voronoi_diagram_traits> VD;
typedef boost::polygon::voronoi_diagram<double> VD;


namespace polytope
{

template<typename RealType>
class BoostTessellator: public Tessellator<2, RealType>
{
public:

  // Constructor, destructor.
  BoostTessellator();
  ~BoostTessellator();

  // Tessellate the given generators. A bounding box is constructed about
  // the generators, and the corners of the bounding box are added as 
  // additional generators if they are not present in the list.
  void tessellate(const std::vector<RealType>& points,
                  Tessellation<2, RealType>& mesh) const;

  // Tessellate with a bounding box representing the boundaries.
  void tessellate(const std::vector<RealType>& points,
                  RealType* low,
                  RealType* high,
                  Tessellation<2, RealType>& mesh) const;

  // Tessellate obeying the given boundaries.
  void tessellate(const std::vector<RealType>& points,
                  const std::vector<RealType>& PLCpoints,
                  const PLC<2, RealType>& geometry,
                  Tessellation<2, RealType>& mesh) const;

  // This Tessellator handles PLCs!
  bool handlesPLCs() const { return true; }

  // The name of the tessellator
  std::string name() const { return "BoostTessellator"; }

  //! Returns the accuracy to which this tessellator can distinguish coordinates.
  //! Should be returned appropriately for normalized coordinates, i.e., if all
  //! coordinates are in the range xi \in [0,1], what is the minimum allowed 
  //! delta in x.
  virtual RealType degeneracy() const { return 1.0e-8; }

private:
  //-------------------- Private interface ---------------------- //

  // ------------------------------------------------- //
  // Specialized tessellations based on the point set  //
  // ------------------------------------------------- //

  // Compute the nodes around a collection of generators
  void computeCellNodes(const std::vector<RealType>& points,
                        std::map<PointType, std::pair<int, int> >& nodeMap,
                        std::vector<std::vector<unsigned> >& cellNodes) const;

  // Compute bounded cell rings from Boost Voronoi diagram
  void computeCellRings(const std::vector<RealType>& points,
			const std::vector<RealType>& PLCpoints,
			const PLC<2, RealType>& geometry,
                        const std::map<PointType, std::pair<int, int> >& nodeMap,
                        std::vector<std::vector<unsigned> >& cellNodes,
			Clipper2d<RealType>& clipper,
                        std::vector<IntRing>& cellRings,
                        const bool collinear,
			const bool performCellAdoption) const;

  // Compute an unbounded tessellation
  void computeVoronoiUnbounded(const std::vector<RealType>& points,
			       Tessellation<2, RealType>& mesh) const;

  // Compute a bounded tessellation
  void computeVoronoiBounded(const std::vector<RealType>& points,
			     const std::vector<RealType>& PLCpoints,
			     const PLC<2, RealType>& geometry,
			     Tessellation<2, RealType>& mesh) const;

  // ----------------------------------------------------- //
  // Private tessellate calls used by internal algorithms  //
  // ----------------------------------------------------- //

  // Bounded tessellation with prescribed bounding box
  void tessellate(const std::vector<RealType>& points,
                  const std::vector<CoordHash>& IntPLCpoints,
                  const PLC<2, RealType>& geometry,
                  const QuantizedCoordinates<2,RealType>& coords,
                  std::vector<std::vector<std::vector<CoordHash> > >& IntCells) const;

  // -------------------------- //
  // Private member variables   //
  // -------------------------- //

  // The quantized coordinates for this tessellator (inner and outer)
  mutable QuantizedCoordinates<2,RealType> mCoords;

  friend class BoostOrphanage<RealType>;

};

} //end polytope namespace

#endif
#endif
