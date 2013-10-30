#ifndef POLYTOPE_POLYGON_HH
#define POLYTOPE_POLYGON_HH

#include "polytope.hh"

namespace polytope
{

template <typename RealType>
class Polygon
{
  public:

  //! Construct a polygon by traversing the given vertices in order.
  //! Here, we assume that vertices[2*i] and vertices[2*i+1] are the 
  //! x and y coordinates of the ith vertex.
  explicit Polygon(const std::vector<RealType>& vertices):
    m_vertices(vertices) 
  {
    POLY_ASSERT((vertices.size() % 2) == 0);
    POLY_ASSERT((vertices.size()/2) >= 3); // Must be at least a triangle!
  }

  //! Construct a polygon by traversing the given vertices in order.
  //! Here, we assume that begin and end are iterators that traverse a 
  //! sequence of x, y coordinates for each vertex in vertex-major order.
  template <typename IteratorType>
  Polygon(IteratorType begin, IteratorType end):
    m_vertices(begin, end) 
  {
    POLY_ASSERT((m_vertices.size() % 2) == 0);
    POLY_ASSERT((m_vertices.size()/2) >= 3); // Must be at least a triangle!
  }

  // Default constructor, copy constructor, destructor, 
  // and assignment operator are generated by the compiler.

  //! Returns the number of vertices in the polygon.
  int numVertices() const { return m_vertices.size()/2; }

  //! Returns the number of edges in the polygon.
  int numEdges() const { return m_vertices.size()/2; }

  //! This enum provides information about the location of a point 
  //! with respect to the PLC.
  enum PointLocation { POINT_INSIDE, POINT_ON_VERTEX, 
                       POINT_ON_EDGE, POINT_OUTSIDE };

  //! Query the given point to determine whether it is inside, outside, 
  //! or on an edge or vertex of this Polygon. This query uses a ray crossing
  //! technique from O'Rourke's "Computational Geometry in C" book.
  PointLocation queryPoint(const RealType point[]) const
  {
    // If the polygon is empty, every point is outside.
    if (m_vertices.empty()) 
      return POINT_OUTSIDE;

    // Numbers of right and left crossings.
    int Rcross = 0, Lcross = 0;

    // Loop over vertices.
    size_t N = m_vertices.size()/2;
    for (size_t i = 0; i < N; ++i)
    {
      // We consider the vertex's coordinates recentered 
      // with the point as the origin.
      RealType vx = m_vertices[2*i]   - point[0],
               vy = m_vertices[2*i+1] - point[1];

      // Check to see whether this point falls on a vertex.
      if ((vx == 0.0) and (vy == 0.0))
        return POINT_ON_VERTEX;

      // Check to see whether the edge e = (i,i+1) straddles 
      // the x axis, with bias above/below.
      size_t i1 = (i+N-1) % N;
      RealType v1x = m_vertices[2*i1]   - point[0],
               v1y = m_vertices[2*i1+1] - point[1];
      bool Rstrad = ((vy > 0.0) != (v1y > 0.0));
      bool Lstrad = ((vy < 0.0) != (v1y < 0.0));

      if (Rstrad or Lstrad) 
      {
        // Compute the intersection of e with the x axis. 
        RealType x = (vx * v1y - v1x * vy) / (v1y - vy);

        if (Rstrad and (x > 0)) 
          Rcross++;
        if (Lstrad and (x < 0))
          Lcross++;
      }
    }

    // The point sits on an edge if Rcross and Lcross 
    // have different parities.
    if ((Rcross % 2) != (Lcross % 2))
      return POINT_ON_EDGE;

    // The point is inside the polygon if the crossings are odd.
    if ((Rcross % 2) == 1)
      return POINT_INSIDE;
    else
      return POINT_OUTSIDE;
  }

  //! Finds the nearest edge in the polygon to the given point, or returns
  //! -1 if the polygon is empty (has no vertices).
  int nearestEdge(const RealType point[]) const
  {
    RealType minDist = FLT_MAX;
    int edge = -1;
    for (int e = 0; e < numEdges(); ++e)
    {
      RealType dist = distanceToEdge(point, e);
      if (dist < minDist)
      {
        minDist = dist;
        edge = e;
      }
    }
    return edge;
  }

  //! Return the mid-point of the nearest edge.
  void nearestEdgePosition(const RealType point[],
                           RealType result[]) const
  {
    const int i0 = nearestEdge(point);
    const int i1 = (i0 + 1) % m_vertices.size();
    result[0] = 0.5*(m_vertices[2*i0    ] + m_vertices[2*i1    ]);
    result[1] = 0.5*(m_vertices[2*i0 + 1] + m_vertices[2*i1 + 1]);
  }

  //! Projects the given point to the given edge of the polygon.
  //! Here, the ith edge contains the ith and ((i+1)%N)th 
  //! vertices of the polygon.
  void project(const RealType point[], int edge, RealType projection[]) const
  {
    POLY_ASSERT(edge >= 0);
    POLY_ASSERT(edge < numEdges());

    // Retrieve the coordinates of the edge's vertices.
    RealType v1x = m_vertices[2*edge], v1y = m_vertices[2*edge+1],
             v2x = m_vertices[2*(edge+1)], v2y = m_vertices[2*(edge+1)+1];

    // Compute the square length of the edge.
    RealType l2 = (v2x-v1x)*(v2x-v1x) + (v2y-v2y)*(v2y-v2y);
    if (l2 == 0.0) // Degenerate edge case.
    {
      // The point is projected to the first vertex.
      projection[0] = v1x;
      projection[1] = v1y;
    }

    // Consider the line extending the segment, parameterized as v1 + t (v2 - v1).
    // We find the projection of the point p onto the line. 
    // It falls where t = [(p-v1) o (v2-v1)] / |v2-v1|^2
    RealType t = ((point[0]-v1x)*(v2x-v1x) + (point[1]-v1y)*(v2y-v1y)) / l2;
    if (t < 0.0)
      t = 0.0;
    else if (t > 1.0)
      t = 1.0;

    // Compute the projection.
    projection[0] = v1x + t * (v2x - v1x);
    projection[1] = v1y + t * (v2y - v1y);
  }

  //! Projects the given point to the nearest edge of the polygon.
  void project(const RealType point[], RealType projection[]) const
  {
    RealType minDist = FLT_MAX;
    for (int e = 0; e < numEdges(); ++e)
    {
      RealType proj[2];
      project(point, e, proj);
      RealType dist = std::sqrt((point[0]-proj[0])*(point[0]-proj[0]) + 
                                (point[1]-proj[1])*(point[1]-proj[1]));
      if (dist < minDist)
      {
        minDist = dist;
        projection[0] = proj[0];
        projection[1] = proj[1];
      }
    }
  }

  //! Computes the minimum distance from this point to the given edge 
  //! in the polygon. Here, the ith edge contains the ith and ((i+1)%N)th 
  //! vertices of the polygon.
  RealType distanceToEdge(const RealType point[], const int i0) const
  {
    // Compute the projection of the point to the given edge.
    const int i1 = (i0 + 1) % m_vertices.size();
    RealType proj[2] = {0.5*(m_vertices[2*i0    ] + m_vertices[2*i1    ]),
                        0.5*(m_vertices[2*i0 + 1] + m_vertices[2*i1 + 1])};

    // Return the distance to the projection.
    return std::sqrt((point[0]-proj[0])*(point[0]-proj[0]) + 
                     (point[1]-proj[1])*(point[1]-proj[1]));
  }


  //! Find the intersection of the given line segement (s1->s2) with the polygon,
  //! closest to the line segments first point (s1).
  //! Returns the edge index for the edge intersected (-1 if none).
  int closestIntersection(const RealType* s1, const RealType* s2,
                          RealType* intersect) const
  {
    const double sx = s2[0] - s1[0];
    const double sy = s2[1] - s1[1];
    unsigned i, j, n = m_vertices.size()/2;
    double ex, ey, esx, esy, phia, phib, thpt, minPhi = FLT_MAX;
    int result = -1;
    for (i = 0; i != n; ++i)
    {
      j = (i + 1) % n;
      ex = m_vertices[2*j] - m_vertices[2*i];
      ey = m_vertices[2*j+1] - m_vertices[2*i+1];
      esx = s1[0] - m_vertices[2*i];
      esy = s1[1] - m_vertices[2*i+1];
      thpt = ey*sx - ex*sy;
      if (thpt != 0.0) {
        phia = (ex*esy - ey*esx)/thpt;
        phib = (sx*esy - sy*esx)/thpt;
        if (phia >= 0.0 and phia <= 1.0 and
            phib >= 0.0 and phib <= 1.0 and
            phia < minPhi)
        {
          minPhi = phia;
          result = i;
          intersect[0] = s1[0] + phia*sx;
          intersect[1] = s1[1] + phia*sy;
        }
      }
    }
    return result;
  }

  // Allow const access to the vertices.
  const std::vector<RealType>& vertices() const { return m_vertices; }

  private:

  //! An array containing the vertices of the polygon. The array has 2N 
  //! entries, where N is the number of vertices, and vertices[2*i] and 
  //! vertices[2*i+1] are the x and y coordinates of the ith vertex.
  std::vector<RealType> m_vertices;
};

}

#endif
