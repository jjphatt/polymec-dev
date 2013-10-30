//------------------------------------------------------------------------------
// Thinly adapted from the following CSG libraries:
//
// Original CSG.JS library by Evan Wallace (http://madebyevan.com), under the MIT license.
// GitHub: https://github.com/evanw/csg.js/
// 
// C++ port by Tomasz Dabrowski (http://28byteslater.com), under the MIT license.
// GitHub: https://github.com/dabroz/csgjs-cpp/
// 
// We have adapted the above for use the the polytope::ReducedPLC classes in 3D.
// Perhaps we could similarly produce the 2D version to provide an option to
// Boost.Geometry?
//------------------------------------------------------------------------------
#ifndef __Polytope_PLC_CSG_3d__
#define __Polytope_PLC_CSG_3d__

#include <map>
#include <vector>
#include <algorithm>

#include "ReducedPLC.hh"
#include "Point.hh"
#include "polytope_geometric_utilities.hh"

namespace polytope {
namespace CSG {

//------------------------------------------------------------------------------
// Public interface methods (forward declaration).
//------------------------------------------------------------------------------
template<typename RealType> ReducedPLC<3, RealType> csg_union    (const ReducedPLC<3, RealType>& a, const ReducedPLC<3, RealType>& b);
template<typename RealType> ReducedPLC<3, RealType> csg_intersect(const ReducedPLC<3, RealType>& a, const ReducedPLC<3, RealType>& b);
template<typename RealType> ReducedPLC<3, RealType> csg_subtract (const ReducedPLC<3, RealType>& a, const ReducedPLC<3, RealType>& b);

//------------------------------------------------------------------------------
// Everything from here down is implementation detail.
//
// PointType vector methods.
//------------------------------------------------------------------------------
namespace CSG_internal_3d {
template<typename RealType> Point3<RealType> lerp(const Point3<RealType>& a, const Point3<RealType>& b, const double v) { 
  const Point3<double> ra(a.x, a.y, a.z), 
                       rb(b.x, b.y, b.z),
                       rc = ra + (rb - ra)*v;
  return Point3<RealType>(RealType(rc.x), RealType(rc.y), RealType(rc.z));
}

// template<typename RealType> Point3<RealType> lerp(const Point3<RealType>& a, const Point3<RealType>& b, const double v) { return a + (b - a) * v; }
template<typename RealType> RealType         length2(const Point3<RealType>& a) { return a.x*a.x + a.y*a.y + a.z*a.z; }
// template<typename RealType> Point3<RealType> unit(const Point3<RealType>& a) { Point3<RealType> result(a); geometry::unitVector<3, RealType>(&result.x); return result; }
template<typename RealType> RealType         dot(const Point3<RealType>& a, const Point3<RealType>& b) { return geometry::dot<3, RealType>(&a.x, &b.x); }
template<typename RealType> Point3<RealType> cross(const Point3<RealType>& a, const Point3<RealType>& b) { Point3<RealType> result; geometry::cross<3, RealType>(&a.x, &b.x, &result.x); return result; }

//------------------------------------------------------------------------------
// Represents a vertex of a polygon. 
//------------------------------------------------------------------------------
template<typename RealType>
struct Vertex {
  typedef Point3<RealType> PointType;
  PointType pos;
  Vertex(): pos() {}
  Vertex(const PointType& posi): pos(posi) {}
};

//------------------------------------------------------------------------------
// Forward declarations.
//------------------------------------------------------------------------------
template<typename RealType> struct Plane;
template<typename RealType> struct Polygon;
template<typename RealType> struct Node;

//------------------------------------------------------------------------------
// Represents a plane in 3D space.
//------------------------------------------------------------------------------
template<typename RealType>
struct Plane {
  typedef Point3<RealType> PointType;
  PointType normal;
  double w;
  static double EPSILON;

  Plane(): normal(), w(0.0) {}
  Plane(const PointType& a, const PointType& b, const PointType& c): 
    normal(cross(b - a, c - a)),
    w(dot(normal, a)) {}
  bool ok() const { return length2(this->normal) > 0.0; }
  void flip() { normal = -normal; w *= -1; }
  void splitPolygon(const Polygon<RealType> & polygon,
                    std::vector<Polygon<RealType> > & coplanarFront,
                    std::vector<Polygon<RealType> > & coplanarBack,
                    std::vector<Polygon<RealType> > & front, 
                    std::vector<Polygon<RealType> > & back) const;
};

//------------------------------------------------------------------------------
// Represents a convex polygon. The vertices used to initialize a polygon must
// be coplanar and form a convex loop.
//------------------------------------------------------------------------------
template<typename RealType>
struct Polygon {
  std::vector<Vertex<RealType> > vertices;
  Plane<RealType> plane;

  void flip() {
    std::reverse(vertices.begin(), vertices.end());
    plane.flip();
  }

  Polygon() {};
  Polygon(const std::vector<Vertex<RealType> > & list):
    vertices(list), 
    plane(vertices[0].pos, vertices[1].pos, vertices[2].pos) {}
};

//------------------------------------------------------------------------------
// Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
// by picking a polygon to split along. That polygon (and all other coplanar
// polygons) are added directly to that node and the other polygons are added to
// the front and/or back subtrees. This is not a leafy BSP tree since there is
// no distinction between internal and leaf nodes.
//------------------------------------------------------------------------------
template<typename RealType> 
struct Node {
  std::vector<Polygon<RealType> > polygons;
  Node<RealType> * front;
  Node<RealType> * back;
  Plane<RealType> plane;

  Node(): front(0), back(0) {}
  Node(const std::vector<Polygon<RealType> > & list);
  ~Node() { delete front; delete back; }

  Node * clone() const;
  void clipTo(const Node * other);
  void invert();
  void build(const std::vector<Polygon<RealType> > & polygon);
  std::vector<Polygon<RealType> > clipPolygons(const std::vector<Polygon<RealType> > & list) const;
  std::vector<Polygon<RealType> > allPolygons() const;
};

//------------------------------------------------------------------------------
// Vertex implementation
//------------------------------------------------------------------------------
// Invert all orientation-specific data. Called when the orientation of a polygon is flipped.
template<typename RealType> Vertex<RealType> flip(const Vertex<RealType>& v) { return v; }

// Create a new vertex between this vertex and `other` by linearly
// interpolating all properties using a parameter of `t`. Subclasses should
// override this to interpolate additional properties.
template<typename RealType>
Vertex<RealType> interpolate(const Vertex<RealType>& a, const Vertex<RealType>& b, double t) {
  Vertex<RealType> ret;
  ret.pos = lerp(a.pos, b.pos, t);
  return ret;
}

//------------------------------------------------------------------------------
// Plane implementation
//------------------------------------------------------------------------------
// Split `polygon` by this plane if needed, then put the polygon or polygon
// fragments in the appropriate lists. Coplanar polygons go into either
// `coplanarFront` or `coplanarBack` depending on their orientation with
// respect to this plane. Polygons in front or in back of this plane go into
// either `front` or `back`.
template<typename RealType>
void Plane<RealType>::splitPolygon(const Polygon<RealType>& polygon, 
                                   std::vector<Polygon<RealType> >& coplanarFront, 
                                   std::vector<Polygon<RealType> >& coplanarBack,
                                   std::vector<Polygon<RealType> > & front, 
                                   std::vector<Polygon<RealType> > & back) const {
  enum {
    COPLANAR = 0,
    FRONT = 1,
    BACK = 2,
    SPANNING = 3
  };

  // Classify each point as well as the entire polygon into one of the above
  // four classes.
  int polygonType = 0;
  std::vector<int> types;

  for (size_t i = 0; i < polygon.vertices.size(); i++) {
    double t = dot(this->normal, polygon.vertices[i].pos) - this->w;
    int type = (t < -EPSILON) ? BACK : ((t > EPSILON) ? FRONT : COPLANAR);
    polygonType |= type;
    types.push_back(type);
  }

  // Put the polygon in the correct list, splitting it when necessary.
  switch (polygonType) {
  case COPLANAR: 
    {
      if (dot(this->normal, polygon.plane.normal) > 0)
        coplanarFront.push_back(polygon);
      else 
        coplanarBack.push_back(polygon);
      break;
    }
  case FRONT: 
    {
      front.push_back(polygon);
      break;
    }
  case BACK: 
    {
      back.push_back(polygon);
      break;
    }
  case SPANNING: 
    {
      std::vector<Vertex<RealType> > f, b;
      for (size_t i = 0; i < polygon.vertices.size(); i++) {
        int j = (i + 1) % polygon.vertices.size();
        int ti = types[i], tj = types[j];
        Vertex<RealType> vi = polygon.vertices[i], vj = polygon.vertices[j];
        if (ti != BACK) f.push_back(vi);
        if (ti != FRONT) b.push_back(vi);
        if ((ti | tj) == SPANNING) {
          double t = (this->w - dot(this->normal, vi.pos)) / dot(this->normal, vj.pos - vi.pos);
          Vertex<RealType> v = interpolate(vi, vj, t);
          f.push_back(v);
          b.push_back(v);
        }
      }
      if (f.size() >= 3) front.push_back(Polygon<RealType>(f));
      if (b.size() >= 3) back.push_back(Polygon<RealType>(b));
      break;
    }
  }
}

//------------------------------------------------------------------------------
// Node implementation
//------------------------------------------------------------------------------
// Return a new CSG solid representing space in either this solid or in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
template<typename RealType>
Node<RealType>* csg_union(const Node<RealType>* a1, const Node<RealType>* b1) {
  Node<RealType> * a = a1->clone();
  Node<RealType> * b = b1->clone();
  a->clipTo(b);
  b->clipTo(a);
  b->invert();
  b->clipTo(a);
  b->invert();
  a->build(b->allPolygons());
  Node<RealType>* ret = new Node<RealType>(a->allPolygons());
  delete a; a = 0;
  delete b; b = 0;
  return ret;
}

// Return a new CSG solid representing space in this solid but not in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
template<typename RealType>
Node<RealType>* csg_subtract(const Node<RealType>* a1, const Node<RealType> * b1) {
  Node<RealType>* a = a1->clone();
  Node<RealType>* b = b1->clone();
  a->invert();
  a->clipTo(b);
  b->clipTo(a);
  b->invert();
  b->clipTo(a);
  b->invert();
  a->build(b->allPolygons());
  a->invert();
  Node<RealType>* ret = new Node<RealType>(a->allPolygons());
  delete a; a = 0;
  delete b; b = 0;
  return ret;
}

// Return a new CSG solid representing space both this solid and in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
template<typename RealType>
Node<RealType>* csg_intersect(const Node<RealType>* a1, const Node<RealType>* b1) {
  Node<RealType>* a = a1->clone();
  Node<RealType>* b = b1->clone();
  a->invert();
  b->clipTo(a);
  b->invert();
  a->clipTo(b);
  b->clipTo(a);
  a->build(b->allPolygons());
  a->invert();
  Node<RealType>* ret = new Node<RealType>(a->allPolygons());
  delete a; a = 0;
  delete b; b = 0;
  return ret;
}

// Convert solid space to empty space and empty space to solid space.
template<typename RealType>
void Node<RealType>::invert() {	
  for (size_t i = 0; i < this->polygons.size(); i++) this->polygons[i].flip();
  this->plane.flip();
  if (this->front) this->front->invert();
  if (this->back) this->back->invert();
  std::swap(this->front, this->back);
}

// Recursively remove all polygons in `polygons` that are inside this BSP
// tree.
template<typename RealType>
std::vector<Polygon<RealType> > 
Node<RealType>::clipPolygons(const std::vector<Polygon<RealType> >& list) const {
  if (!this->plane.ok()) return list;
  std::vector<Polygon<RealType> > list_front, list_back;	
  for (size_t i = 0; i < list.size(); i++) {
    this->plane.splitPolygon(list[i], list_front, list_back, list_front, list_back);
  }
  if (this->front) list_front = this->front->clipPolygons(list_front);
  if (this->back) list_back = this->back->clipPolygons(list_back);
  else list_back.clear();
	
  list_front.insert(list_front.end(), list_back.begin(), list_back.end());
  return list_front;
}

// Remove all polygons in this BSP tree that are inside the other BSP tree
// `bsp`.
template<typename RealType>
void Node<RealType>::clipTo(const Node<RealType>* other) {
  this->polygons = other->clipPolygons(this->polygons);
  if (this->front) this->front->clipTo(other);
  if (this->back) this->back->clipTo(other);
}

// Return a list of all polygons in this BSP tree.
template<typename RealType>
std::vector<Polygon<RealType> > Node<RealType>::allPolygons() const {
  std::vector<Polygon<RealType> > list = this->polygons;
  std::vector<Polygon<RealType> > list_front, list_back;
  if (this->front) list_front = this->front->allPolygons();
  if (this->back) list_back = this->back->allPolygons();
  list.insert(list.end(), list_front.begin(), list_front.end());
  list.insert(list.end(), list_back.begin(), list_back.end());
  return list;
}

template<typename RealType>
Node<RealType>* Node<RealType>::clone() const {
  Node<RealType>* ret = new Node<RealType>();
  ret->polygons = this->polygons;
  ret->plane = this->plane;
  if (this->front) ret->front = this->front->clone();
  if (this->back) ret->back = this->back->clone();
  return ret;
}

// Build a BSP tree out of `polygons`. When called on an existing tree, the
// new polygons are filtered down to the bottom of the tree and become new
// nodes there. Each set of polygons is partitioned using the first polygon
// (no heuristic is used to pick a good split).
template<typename RealType>
void Node<RealType>::build(const std::vector<Polygon<RealType> > & list) {
  if (!list.size()) return;
  if (!this->plane.ok()) this->plane = list[0].plane;
  std::vector<Polygon<RealType> > list_front, list_back;
  for (size_t i = 0; i < list.size(); i++) {
    this->plane.splitPolygon(list[i], this->polygons, this->polygons, list_front, list_back);
  }
  if (list_front.size()) {
    if (!this->front) this->front = new Node<RealType>;
    this->front->build(list_front);
  }
  if (list_back.size()) {
    if (!this->back) this->back = new Node<RealType>;
    this->back->build(list_back);
  }
}

template<typename RealType>
Node<RealType>::Node(const std::vector<Polygon<RealType> > & list): front(0), back(0) {
  build(list);
}

//------------------------------------------------------------------------------
// Convert between polytope::ReducedPLC <-> vector<Polygon>
//------------------------------------------------------------------------------
// Convert a ReducedPLC to a set of CSG polygons.
// Note PLC facets are arbitrary polygons, but this CSG code is built with 
// triangular facets.
template<typename RealType>
std::vector<Polygon<RealType> >
ReducedPLCtoPolygons(const ReducedPLC<3, RealType>& model) {
  typedef Point3<RealType> PointType;
  const std::vector<RealType>& coords = model.points;
  std::vector<Polygon<RealType> > list;
  std::vector<Vertex<RealType> > triangle(3);
  const unsigned nfacets = model.facets.size();
  POLY_ASSERT(nfacets >= 4);
  for (size_t i = 0; i != nfacets; ++i) {
    const unsigned n = model.facets[i].size();
    POLY_ASSERT(n >= 3);
    const unsigned k1 = model.facets[i][0];
    for (size_t j = 1; j != n-1; ++j) {
      const unsigned k2 = model.facets[i][j],
                     k3 = model.facets[i][(j + 1) % n];
      triangle[0] = Vertex<RealType>(PointType(coords[3*k1], coords[3*k1 + 1], coords[3*k1 + 2]));
      triangle[1] = Vertex<RealType>(PointType(coords[3*k2], coords[3*k2 + 1], coords[3*k2 + 2]));
      triangle[2] = Vertex<RealType>(PointType(coords[3*k3], coords[3*k3 + 1], coords[3*k3 + 2]));
      list.push_back(Polygon<RealType>(triangle));
    }
  }
  return list;
}

// Convert a set of polygons to a ReducedPLC.
// Note we do not remove degeneracies here, so it's up to the caller
// to do with that as they will!
template<typename RealType>
ReducedPLC<3, RealType>
ReducedPLCfromPolygons(const std::vector<Polygon<RealType> >& polygons) {
  typedef Point3<RealType> PointType;
  ReducedPLC<3, RealType> result;
  const unsigned nfacets = polygons.size();
  POLY_ASSERT(nfacets >= 4);
  for (size_t i = 0; i != nfacets; ++i) {
    // std::cerr << "Polygon " << i << " : " << polygons[i].vertices[0].pos << " " << polygons[i].vertices[1].pos;
    const unsigned n = polygons[i].vertices.size();
    for (size_t j = 2; j < n; ++j) {
      // std::cerr << " " << polygons[i].vertices[j].pos;
      result.facets.push_back(std::vector<int>());
      result.facets.back().push_back(result.points.size()/3);
      result.points.push_back(polygons[i].vertices[0].pos.x);
      result.points.push_back(polygons[i].vertices[0].pos.y);
      result.points.push_back(polygons[i].vertices[0].pos.z);
      result.facets.back().push_back(result.points.size()/3);
      result.points.push_back(polygons[i].vertices[j - 1].pos.x);
      result.points.push_back(polygons[i].vertices[j - 1].pos.y);
      result.points.push_back(polygons[i].vertices[j - 1].pos.z);
      result.facets.back().push_back(result.points.size()/3);
      result.points.push_back(polygons[i].vertices[j].pos.x);
      result.points.push_back(polygons[i].vertices[j].pos.y);
      result.points.push_back(polygons[i].vertices[j].pos.z);
    }
    // std::cerr << std::endl;
  }
  return result;
}

} // end namespace CSG_internal

//------------------------------------------------------------------------------
// Public interface implementation
//------------------------------------------------------------------------------
template<typename RealType>
ReducedPLC<3, RealType> csg_union(const ReducedPLC<3, RealType>& a,
                                  const ReducedPLC<3, RealType>& b) {
  using namespace CSG_internal_3d;
  Node<RealType>* A = new Node<RealType>(ReducedPLCtoPolygons(a));
  Node<RealType>* B = new Node<RealType>(ReducedPLCtoPolygons(b));
  Node<RealType>* AB = csg_union(A, B);                                  // <-- only difference
  std::vector<Polygon<RealType> > polygons = AB->allPolygons();
  delete A; A = 0;
  delete B; B = 0;
  delete AB; AB = 0;
  return ReducedPLCfromPolygons(polygons);
}

template<typename RealType>
ReducedPLC<3, RealType> csg_intersect(const ReducedPLC<3, RealType>& a,
                                      const ReducedPLC<3, RealType>& b) {
  using namespace CSG_internal_3d;
  Node<RealType>* A = new Node<RealType>(ReducedPLCtoPolygons(a));
  Node<RealType>* B = new Node<RealType>(ReducedPLCtoPolygons(b));
  Node<RealType>* AB = csg_intersect(A, B);                           // <-- only difference
  std::vector<Polygon<RealType> > polygons = AB->allPolygons();
  delete A; A = 0;
  delete B; B = 0;
  delete AB; AB = 0;
  return ReducedPLCfromPolygons(polygons);
}

template<typename RealType>
ReducedPLC<3, RealType> csg_subtract(const ReducedPLC<3, RealType>& a,
                                     const ReducedPLC<3, RealType>& b) {
  using namespace CSG_internal_3d;
  Node<RealType>* A = new Node<RealType>(ReducedPLCtoPolygons(a));
  Node<RealType>* B = new Node<RealType>(ReducedPLCtoPolygons(b));
  Node<RealType>* AB = csg_subtract(A, B);                            // <-- only difference
  std::vector<Polygon<RealType> > polygons = AB->allPolygons();
  delete A; A = 0;
  delete B; B = 0;
  delete AB; AB = 0;
  return ReducedPLCfromPolygons(polygons);
}

} // end namespace CGS
} // end namespace polytope

#endif
