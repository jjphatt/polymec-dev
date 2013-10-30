// TriangleUnbounded
// Try to tessellate a series of degenerate unbounded tessellations and
// check that the correct number of cells, nodes, and faces result

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

#define POLY_CHECK_BOOL(x) if (!(x)) { return false; }

using namespace std;
using namespace polytope;


// -----------------------------------------------------------------------
// checkMesh
// -----------------------------------------------------------------------
bool checkMesh(const Tessellation<2,double>& mesh,
               const unsigned ncells,
               const unsigned nnodes,
               const unsigned nfaces,
               const unsigned ninfNodes,
               const unsigned ninfFaces) {
  POLY_CHECK_BOOL(mesh.cells.size()   == ncells);
  POLY_CHECK_BOOL(mesh.nodes.size()/2 == nnodes);
  POLY_CHECK_BOOL(mesh.faces.size()   == nfaces);
  unsigned infNodeCount=0;
  for(unsigned i = 0; i < mesh.infNodes.size(); ++i){
    if( mesh.infNodes[i] == 1 ) ++infNodeCount;
  }
  POLY_CHECK_BOOL(infNodeCount == ninfNodes);
  unsigned infFaceCount=0;
  for(unsigned i = 0; i < mesh.infFaces.size(); ++i){
    if( mesh.infFaces[i] == 1 ) ++infFaceCount;
  }
  POLY_CHECK_BOOL(infFaceCount == ninfFaces);
  return true;
}


// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator) {

  string testName = "Unbounded_" + tessellator.name();
  int ntest = 1;
  vector<int> howDidIDo;
  Tessellation<2,double> mesh;

  // Test 1: Circle of generators
  {
    cout << "\nTest " << ntest << ": Circle of generators" << endl;
    int N = 18;
    vector<double> points(2*N);
    for (unsigned i = 0; i < N; ++i){
       double theta = 2.0*M_PI*double(i)/double(N);
       points[2*i] = cos(theta);  points[2*i+1] = sin(theta);
    }
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, ntest);
    bool pass = checkMesh(mesh,N,N+1,2*N,N,N);
    if (!pass) howDidIDo.push_back(ntest);
    ++ntest;
  }

  // Test 2: Circle of generators, random center
  {
    cout << "\nTest " << ntest << ": Circle of generators, random center" << endl;
    int N = 18;
    vector<double> points(2*N);
    double center[2] = {random01(), random01()};
    for (unsigned i = 0; i < N; ++i){
       double theta = 2.0*M_PI*double(i)/double(N);
       points[2*i  ] = center[0] + cos(theta);
       points[2*i+1] = center[1] + sin(theta);
    }
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, ntest);
    bool pass = checkMesh(mesh,N,N+1,2*N,N,N);
    if (!pass) howDidIDo.push_back(ntest);
    ++ntest;
  }

  // Test 3: Two uniform rows of generators
  {
    cout << "\nTest " << ntest << ": Two uniform rows of generators" << endl;
    int N = 10;
    vector<double> points(4*N);
    for (unsigned i = 0; i < N; ++i){
       points[4*i  ] = i;  points[4*i+1] = -1.0;
       points[4*i+2] = i;  points[4*i+3] =  1.0;
    }
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, ntest);
    bool pass = checkMesh(mesh, 2*N, 3*N-1, 5*N-2, 2*N, 2*N);
    if (!pass) howDidIDo.push_back(ntest);
    ++ntest;
  }

  // Test 4: Collinear generators with one non-collinear
  {
    cout << "\nTest " << ntest << ": Collinear generators, except one" << endl;
    int N = 10;
    vector<double> points(2*N);
    for (unsigned i = 0; i < N; ++i)  points[2*i] = double(i);
    points.push_back(4.5);
    points.push_back(1.0);
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, ntest);
    ++ntest;
  }

  // Test 5: 2x2 Cartesian Generators
  {
    cout << "\nTest 4: 2x2 Cartesian generators" << endl;
    vector<double> points;
    points.push_back(0.0); points.push_back(0.0);
    points.push_back(1.0); points.push_back(0.0);
    points.push_back(1.0); points.push_back(1.0);
    points.push_back(0.0); points.push_back(1.0);
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, ntest);
    bool pass = checkMesh(mesh, 4, 5, 8, 4, 4);
    if (!pass) howDidIDo.push_back(ntest);
    ++ntest;
  }

  // Test 6: Two generators
  {
    cout << "\nTest " << ntest << ": Two generators" << endl;
    vector<double> points;
    points.push_back(0.0); points.push_back(0.0);
    points.push_back(1.0); points.push_back(0.0);
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, ntest);
    bool pass = checkMesh(mesh, 2, 4, 5, 4, 4);
    if (!pass) howDidIDo.push_back(ntest);
    ++ntest;
  }

  // Test 7: Line of generators, uniform
  {
    cout << "\nTest " << ntest << ": Uniform line of generators" << endl;
    int N=10;
    vector<double> points(2*N, 0);
    for (unsigned i = 0; i < N; ++i)  points[2*i] = double(i);
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, ntest);
    bool pass = checkMesh(mesh, N, 2*N, 3*N-1, 2*N, 2*N);
    if (!pass) howDidIDo.push_back(ntest);
    ++ntest;
  }

  // Test 8: Line of generators, non-uniform
  {
    cout << "\nTest " << ntest << ": Non-uniform line of generators" << endl;
    int N = 10;
    vector<double> points(2*N);
    for (unsigned i = 0; i < N; ++i)  points[2*i] = double(i) + random01() - 0.5;
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, ntest);
    bool pass = checkMesh(mesh, N, 2*N, 3*N-1, 2*N, 2*N);
    if (!pass) howDidIDo.push_back(ntest);
    ++ntest;
  }

  // Test 9: Line of generators, non-uniform, shuffled
  {
    cout << "\nTest " << ntest << ": Non-uniform line of generators, shuffled" << endl;
    int N = 10;
    vector<double> points(2*N);
    int indices[10] = {5,7,1,2,8,0,3,9,4,6};
    for (unsigned i = 0; i < N; ++i)  points[2*i] = double(indices[i]) + random01() - 0.5;
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, mesh);
    outputMesh(mesh, testName, points, ntest);
    bool pass = checkMesh(mesh, N, 2*N, 3*N-1, 2*N, 2*N);    
    if (!pass) howDidIDo.push_back(ntest);
    ++ntest;
  } 

  if (!howDidIDo.empty()) {
    for (unsigned i = 0; i != howDidIDo.size(); ++i)
      cout << "Failed Test " << howDidIDo[i] << endl;
    POLY_ASSERT(false);
  }
}


// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

#if HAVE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator<double> tessellator;
    test(tessellator);  
  }
#endif   


#if HAVE_BOOST_VORONOI
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    test(tessellator);
  }
#endif
  

  cout << "PASS" << endl;

#if HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
