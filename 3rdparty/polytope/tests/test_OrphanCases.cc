// test_OrphanCases
//
// A collection of difficult test cases for the orphaned cell algorithm.
// All tests involve meshing a circular region with a star-shaped hole
// in the middle with only 20 points. Seeding the random number generator
// provides the input generator locations. Test cases were found through 
// trial-and-error (though with obnoxious frequency!).

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "Boundary2D.hh"
#include "Generators.hh"
#include "polytope_test_utilities.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// printArea
// -----------------------------------------------------------------------
void printArea(Boundary2D<double>& boundary,
	       Tessellation<2,double>& mesh) {
   const double area = computeTessellationArea(mesh);
   const double relErr = (boundary.mArea-area)/boundary.mArea;
   const double tol = 1.0e-5;
   cout << "Tessellation Area = " << area << endl;
   cout << "Relative error    = " << std::abs(boundary.mArea-area)/boundary.mArea << endl;
   POLY_ASSERT(relErr < tol);
}

// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator) {

  // Test name for output
  string testName = "OrphanCases_" + tessellator.name();
  
  // Initialize boundary and tessellator
  Boundary2D<double> boundary;
  
  // Circular region with star-shaped hole
  int bType = 5;
  boundary.setDefaultBoundary(bType);
  
  int i = 1;
  
  // Test 1: Cell parents multiple orphans
  {
    srand(10489593);
    cout << "\nTest 1: Cell parents multiple orphans" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(20);
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, generators.mPoints, i);
    printArea(boundary,mesh);
    ++i;
  }
  
  // Test 2: Orphan neighbors are also parents of orphans
  {
    srand(10489594);
    cout << "\nTest 2: Orphan neighbors are also parents of orphans" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(20);
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, generators.mPoints, i);
    printArea(boundary,mesh);
    ++i;
  }
  
  // Test 3: Overlapping orphans
  {
    srand(10489609);
    cout << "\nTest 3: Overlapping orphans" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(20);
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, generators.mPoints, i);
    printArea(boundary,mesh);
    ++i;
  }
  
  // Test 4: Empty orphan neighbor set
  {
    srand(10489611);
    cout << "\nTest 4: Empty orphan neighbor set" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(20);
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, generators.mPoints, i);
    printArea(boundary,mesh);
    ++i;
  }
  
  // Test 5: Boost.Geometry calls invalid overlay exception
  {
    srand(10489612);
    cout << "\nTest 5: Boost.Geometry calls invalid overlay exception" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(20);
    Tessellation<2,double> mesh;
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, generators.mPoints, i);
    printArea(boundary,mesh);
    ++i;
  }
  
  // Test 6: Nonconvex boundary with three internal generators
  {
    cout << "\nTest 6: Nonconvex boundary with three internal generators" << endl;
    std::vector<double> points;
    points.push_back(0.05); points.push_back(0.60);
    points.push_back(0.40); points.push_back(0.10);
    points.push_back(0.60); points.push_back(0.50);
    std::vector<double> PLCpoints;
    PLCpoints.push_back(0.0); PLCpoints.push_back(0.0);
    PLCpoints.push_back(0.1); PLCpoints.push_back(0.0);
    PLCpoints.push_back(0.2); PLCpoints.push_back(0.8);
    PLCpoints.push_back(0.3); PLCpoints.push_back(0.0);
    PLCpoints.push_back(1.0); PLCpoints.push_back(0.0);
    PLCpoints.push_back(1.0); PLCpoints.push_back(1.0);
    PLCpoints.push_back(0.0); PLCpoints.push_back(1.0);
    PLC<2,double> boundary;
    boundary.facets.resize(7, std::vector<int>(2));
    for (int j = 0; j < 7; ++j){
      boundary.facets[j][0] = j;
      boundary.facets[j][1] = (j+1) % 7;
    }
    Tessellation<2,double> mesh;
    tessellator.tessellate(points, PLCpoints, boundary, mesh);
    outputMesh(mesh, testName, points, i);
    ++i;
  }  

  // Test 7: Original 3x3 Test Case
  {
    cout << "\nTest 7: 3x3 Unit Test with 2 Orphans" << endl;
    vector<double> PLCpoints;
    vector<double> points;
    PLC<2,double> boundary;
    Tessellation<2,double> mesh;
    PLCpoints.push_back(0.0);  PLCpoints.push_back(0.0);
    PLCpoints.push_back(1.2);  PLCpoints.push_back(0.0);
    PLCpoints.push_back(1.2);  PLCpoints.push_back(1.3);
    PLCpoints.push_back(1.3);  PLCpoints.push_back(1.3);
    PLCpoints.push_back(1.3);  PLCpoints.push_back(0.0);
    PLCpoints.push_back(3.0);  PLCpoints.push_back(0.0);
    PLCpoints.push_back(3.0);  PLCpoints.push_back(3.0);
    PLCpoints.push_back(1.3);  PLCpoints.push_back(3.0);
    PLCpoints.push_back(1.3);  PLCpoints.push_back(1.7);
    PLCpoints.push_back(1.2);  PLCpoints.push_back(1.7);
    PLCpoints.push_back(1.2);  PLCpoints.push_back(3.0);
    PLCpoints.push_back(0.0);  PLCpoints.push_back(3.0);
    
    int ix, iy, nx = 3;
    double xi, yi;
    for (iy = 0; iy != nx; ++iy) {
      yi = iy + 0.5;
      for (ix = 0; ix != nx; ++ix) {
	xi = ix + 0.5;
	points.push_back(xi);  points.push_back(yi);
      }
    }
    
    int nSides = PLCpoints.size()/2;
    boundary.facets.resize( nSides, std::vector<int>(2) );
    for (unsigned j = 0; j != nSides; ++j){
      boundary.facets[j][0] = j;
      boundary.facets[j][1] = (j+1) % nSides;
    }
    
    tessellator.tessellate(points, PLCpoints, boundary, mesh);
    outputMesh(mesh, testName, points, i);
    const double trueArea = 8.74;
    const double tessArea = computeTessellationArea(mesh);
    const double fracerr  = std::abs(trueArea - tessArea)/trueArea;
    const double tol      = 1.0e-7;
    POLY_CHECK2(fracerr < tol, "Relative error in the tessellation "
        	<< "area exceeds tolerance:" << endl
        	<< "      Area = " << tessArea << endl
        	<< "     Error = " << trueArea - tessArea << endl
        	<< "Frac Error = " << fracerr);
    ++i;
  }

  // Test 8: Lots of random points
  {
    cout << "\nTest 8: Lots of random points" << endl;
    const unsigned N = 100;
    for (unsigned iter = 0; iter != N; ++iter) {
      srand(1049520+iter);
      Generators<2,double> generators(boundary);
      generators.randomPoints(50);
      Tessellation<2,double> mesh;
      tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
      outputMesh(mesh, testName, generators.mPoints, i+iter);
      cout << iter << endl;
      printArea(boundary,mesh);
    }
    ++i;
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
