//----------------------------------------------------------------------------//
// DistributedTessellator
//----------------------------------------------------------------------------//
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>
#include <list>
#include <limits>
#include "mpi.h"

#include "polytope.hh"
#include "Point.hh"
#include "ReducedPLC.hh"
#include "convexHull_2d.hh"
#include "convexHull_3d.hh"
#include "convexIntersect.hh"
#include "deleteCells.hh"
#include "bisectSearch.hh"
#include "polytope_serialize.hh"
#include "polytope_parallel_utilities.hh"
#include "DimensionTraits.hh"
#include "mortonOrderIndices.hh"
#include "checkDistributedTessellation.hh"

using namespace std;
using std::min;
using std::max;
using std::abs;

namespace {  // We hide internal functions in an anonymous namespace

//------------------------------------------------------------------------------
// Comparator to compare std::pair's by their first or second element.
//------------------------------------------------------------------------------
template<typename T1, typename T2>
struct ComparePairByFirstElement {
  bool operator()(const std::pair<T1, T2>& lhs, const std::pair<T1, T2>& rhs) const {
    return lhs.first < rhs.first;
  }
};

template<typename T1, typename T2>
struct ComparePairBySecondElement {
  bool operator()(const std::pair<T1, T2>& lhs, const std::pair<T1, T2>& rhs) const {
    return lhs.second < rhs.second;
  }
};

//------------------------------------------------------------------------------
// Sort a vector of stuff by the given keys.
//------------------------------------------------------------------------------
template<typename Value, typename Key>
void
sortByKeys(vector<Value>& values, const vector<Key>& keys) {
  POLY_ASSERT(values.size() == keys.size());
  vector<pair<Key, Value> > stuff;
  stuff.reserve(values.size());
  for (unsigned i = 0; i != values.size(); ++i) stuff.push_back(make_pair(keys[i], values[i]));
  POLY_ASSERT(stuff.size() == values.size());
  sort(stuff.begin(), stuff.end(), ComparePairByFirstElement<Key, Value>());
  for (unsigned i = 0; i != values.size(); ++i) values[i] = stuff[i].second;
}

} // end anonymous namespace

namespace polytope {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
DistributedTessellator<Dimension, RealType>::
DistributedTessellator(Tessellator<Dimension, RealType>* tessellator,
                       bool assumeControl,
                       bool buildCommunicationInfo):
  mSerialTessellator(tessellator),
  mAssumeControl(assumeControl),
  mBuildCommunicationInfo(buildCommunicationInfo),
  mType(unbounded),
  mLow(0),
  mHigh(0),
  mPLCpointsPtr(0),
  mPLCptr(0) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
DistributedTessellator<Dimension, RealType>::
~DistributedTessellator() {
  if (mAssumeControl)
    delete mSerialTessellator;
}

//------------------------------------------------------------------------------
// Compute an unbounded tessellation.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
tessellate(const vector<RealType>& points,
           Tessellation<Dimension, RealType>& mesh) const {
  mType = unbounded;
  this->computeDistributedTessellation(points, mesh);
}

//------------------------------------------------------------------------------
// Compute the tessellation in the box.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
tessellate(const vector<RealType>& points,
           RealType* low,
           RealType* high,
           Tessellation<Dimension, RealType>& mesh) const {
  mType = box;
  mLow = low;
  mHigh = high;
  this->computeDistributedTessellation(points, mesh);
}

//------------------------------------------------------------------------------
// Compute the tessellation in a PLC.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
tessellate(const vector<RealType>& points,
           const vector<RealType>& PLCpoints,
           const PLC<Dimension, RealType>& geometry,
           Tessellation<Dimension, RealType>& mesh) const {
  mType = plc;
  mPLCpointsPtr = &PLCpoints;
  mPLCptr = &geometry;
  this->computeDistributedTessellation(points, mesh);
}

//------------------------------------------------------------------------------
// This method encapsulates the actual distributed tessellation algorithm.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
computeDistributedTessellation(const vector<RealType>& points,
                               Tessellation<Dimension, RealType>& mesh) const {

  // Some spiffy shorthand typedefs.
  typedef typename DimensionTraits<Dimension, RealType>::ConvexHull ConvexHull;
  typedef typename DimensionTraits<Dimension, RealType>::CoordHash CoordHash;
  typedef typename DimensionTraits<Dimension, RealType>::Point Point;
  typedef typename DimensionTraits<Dimension, RealType>::RealPoint RealPoint;
  typedef KeyTraits::Key Key;
  const double degeneracy = 1.5e-8;
#if HAVE_SILO
  const bool siloBlago = false;
#endif  

  // Parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  const unsigned nlocal = points.size() / Dimension;

  // Compute the bounding box for normalizing our coordinates.
  RealType rlow[Dimension], rhigh[Dimension]; //, genLow[Dimension];
  this->computeBoundingBox(points, rlow, rhigh);
  
  // Copy points to generator vector
  vector<RealType> generators;
  copy(points.begin(), points.end(), back_inserter(generators));
  
  // RealType fscale = 0;
  // for (unsigned i = 0; i != Dimension; ++i) {
  //   fscale = max(fscale, rhigh[i] - rlow[i]);
  //   genLow[i] = RealType(0.0);
  // }
  // POLY_ASSERT(fscale > 0);
  // fscale = 1.0/fscale;

  // // Our goal is to build up the necessary set of generators to completely 
  // // specify the tessellation for all the points passed on this processor.
  // // Start by copying the input set in normalized coordinates.
  // vector<RealType> generators;
  // generators.reserve(points.size());
  // for (unsigned i = 0; i != nlocal; ++i) {
  //   for (unsigned j = 0; j != Dimension; ++j) {
  //     generators.push_back((points[Dimension*i + j] - rlow[j])*fscale);
  //   }
  // }

  // POLY_ASSERT(generators.size() == points.size());
  // POLY_ASSERT2(nlocal == 0 or *min_element(generators.begin(), generators.end()) >= 0,
  //              "Min element out of range: " << *min_element(generators.begin(), generators.end()));
  // POLY_ASSERT2(nlocal == 0 or *max_element(generators.begin(), generators.end()) <= 1,
  //              "Max element out of range: " << *max_element(generators.begin(), generators.end()));

  // We can skip a lot of work if there's only one domain!
  vector<unsigned> gen2domain(nlocal, rank);
  if (numProcs > 1) {

    // Compute the convex hull of each domain and distribute them to all processes.
    vector<ConvexHull> domainHulls; // (numProcs);
    vector<unsigned> domainCellOffset(1, 0);
    {
      ConvexHull localHull = DimensionTraits<Dimension, RealType>::convexHull(generators, rlow, degeneracy);

      // We can skip some work here if the convex hull is of lower dimension than the 
      // problem itself. We determine the set of visible points for a full-dimension
      // hull by computing a local tessellation and collecting all generators that have
      // an exterior face on the local mesh. For a lower dimension hull, every point
      // is visible.
      if (DimensionTraits<Dimension, RealType>::hullDimension(localHull) == Dimension) {

        Tessellation<Dimension, RealType> localMesh;
        mSerialTessellator->tessellate(points, localHull.points, localHull, localMesh);
        
        // Blago!
#if HAVE_SILO
        if (siloBlago)
        {
          vector<double> r2(localMesh.cells.size(), rank);
          vector<double> px(localMesh.cells.size());
          vector<double> py(localMesh.cells.size());
          vector<double> vis(localMesh.cells.size(), 0.0);
          for (unsigned i = 0; i != localMesh.cells.size(); ++i) {
            px[i] = points[2*i];
            py[i] = points[2*i+1];
          }
          for (unsigned i = 0; i != localMesh.faceCells.size(); ++i) {
            if (localMesh.faceCells[i].size() == 1 ) {
              unsigned icell = ((localMesh.faceCells[i][0] < 0) ? 
                                ~localMesh.faceCells[i][0]      : 
                                localMesh.faceCells[i][0]);
              vis[icell] = 1.0;
            }
          }
          map<string, double*> fields, cellFields;
          cellFields["domain"] = &r2[0];
          cellFields["gen_x"] = &px[0];
          cellFields["gen_y"] = &py[0];
          cellFields["visible"] = &vis[0];
          cerr << "Writing local mesh with " << localMesh.cells.size() << endl;
          polytope::SiloWriter<Dimension, RealType>::write(localMesh, fields, fields, 
                                                           fields, cellFields, 
                                                           "test_DistributedTessellator_localMesh");
          MPI_Barrier(MPI_COMM_WORLD);
        }
#endif
        // Blago!

        set<unsigned> exteriorCells;
        for (unsigned i = 0; i != localMesh.faceCells.size(); ++i) {
          if (localMesh.faceCells[i].size() == 1 ) {
            unsigned icell = ((localMesh.faceCells[i][0] < 0) ? 
                              ~localMesh.faceCells[i][0]      : 
                              localMesh.faceCells[i][0]);
            exteriorCells.insert(icell);
          }        
        }
        
        vector<RealPoint> exteriorPoints;
        for (unsigned i = 0; i != localHull.points.size()/Dimension; ++i) {
          exteriorPoints.push_back(DimensionTraits<Dimension, RealType>::constructPoint(&(localHull.points[Dimension*i])));
        }
        
        for (set<unsigned>::const_iterator itr = exteriorCells.begin();
             itr != exteriorCells.end(); ++itr) {
          RealPoint addPoint = DimensionTraits<Dimension, RealType>::constructPoint(&(points[Dimension * (*itr)]));
          typename vector<RealPoint>::iterator it = std::find(exteriorPoints.begin(), exteriorPoints.end(), addPoint);
          if (it == exteriorPoints.end())  exteriorPoints.push_back(addPoint);
        }   
        
        localHull.points.clear();
        for (typename vector<RealPoint>::iterator pointItr = exteriorPoints.begin();
             pointItr != exteriorPoints.end(); ++pointItr) {
           copy(&(*pointItr).x, (&(*pointItr).x + Dimension), back_inserter(localHull.points));
        }
      }
      
      // We have a lower-dimension hull. Every point is visible
      else {
        cerr << "Lower-dimensional data" << endl;
        Tessellation<Dimension, RealType> localMesh;
        mSerialTessellator->tessellate(points, localMesh);
        localHull.points = generators;
      }

      // Serialize and send
      vector<char> localBuffer;
      serialize(localHull, localBuffer);
      for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
        vector<char> buffer = localBuffer;
        unsigned bufSize = localBuffer.size();
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
        buffer.resize(bufSize);
        MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, MPI_COMM_WORLD);
        vector<char>::const_iterator itr = buffer.begin();
        ConvexHull newHull;
        deserialize(newHull, itr, buffer.end());
        POLY_ASSERT(itr == buffer.end());
        domainHulls.push_back(newHull);
        domainCellOffset.push_back(domainCellOffset.back() + domainHulls[sendProc].points.size()/Dimension);
      }
    }
    POLY_ASSERT(domainHulls.size() == numProcs);
    POLY_ASSERT(domainCellOffset.size() == numProcs + 1);

    // // Blago!
    // cerr << "Domain cell offsets : ";
    // copy(domainCellOffset.begin(), domainCellOffset.end(), ostream_iterator<unsigned>(cerr, " "));
    // cerr << endl;
    // cerr << " mLow, mHigh : (" << mLow[0] << " " << mLow[1] << ") (" << mHigh[0] << " " << mHigh[1] << ")" << endl;
    // // Blago!

    // Create a tessellation of the hull vertices for all domains.
    vector<RealType> hullGenerators;
    for (unsigned i = 0; i != numProcs; ++i) {
      copy(domainHulls[i].points.begin(), domainHulls[i].points.end(), back_inserter(hullGenerators));
    }
    POLY_ASSERT(hullGenerators.size()/Dimension == domainCellOffset.back());
    Tessellation<Dimension, RealType> hullMesh;
    this->tessellationWrapper(hullGenerators, hullMesh);

    // PLC<Dimension, RealType> plcTmp;
    // plcTmp.facets = mPLCptr->facets;
    // mSerialTessellator->tessellate(hullGenerators, *mPLCpointsPtr, plcTmp, hullMesh);


    // Blago!
#if HAVE_SILO
    if (siloBlago)
    {
      vector<double> r2(hullMesh.cells.size(), rank);
      vector<double> r3(hullMesh.cells.size());
      vector<double> px(hullMesh.cells.size());
      vector<double> py(hullMesh.cells.size());
      for (unsigned i = 0; i != hullMesh.cells.size(); ++i) {
        unsigned procOwner = bisectSearch(domainCellOffset, i);
        r3[i] = double(procOwner);
        px[i] = hullGenerators[2*i];
        py[i] = hullGenerators[2*i+1];
      }
      map<string, double*> fields, hullCellFields;
      hullCellFields["domain"] = &r2[0];
      hullCellFields["domainOwner"] = &r3[0];
      hullCellFields["gen_x"] = &px[0];
      hullCellFields["gen_y"] = &py[0];
      cerr << "Writing hull mesh with " << hullMesh.cells.size() << endl;
      polytope::SiloWriter<Dimension, RealType>::write(hullMesh, fields, fields, fields, hullCellFields, "test_DistributedTessellator_hullMesh");
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
    // Blago!

    // Find the set of domains we need to communicate with according to two criteria:
    // 1.  Any domain hull that intersects our own.
    // 2.  Any domain hull that has elements adjacent to one of ours in the hullMesh.
    set<unsigned> neighborSet;

    // First any hulls that intersect ours.
    for (unsigned otherProc = 0; otherProc != numProcs; ++otherProc) {
      if (otherProc != rank and
          convexIntersect(domainHulls[otherProc], domainHulls[rank])) neighborSet.insert(otherProc);
    }

    // We need the set of cells that share nodes.
    const vector<set<unsigned> > hullNodeCells = hullMesh.computeNodeCells();
    const vector<set<unsigned> > hullCellToNodes = hullMesh.computeCellToNodes();
    set<unsigned> cellsOfInterest;

    // Now any hulls that have elements adjacent to ours in the hull mesh.
    for (unsigned icell = domainCellOffset[rank]; icell != domainCellOffset[rank + 1]; ++icell) {
      for (set<unsigned>::const_iterator nodeItr1 = hullCellToNodes[icell].begin();
           nodeItr1 != hullCellToNodes[icell].end(); ++nodeItr1){
        for (set<unsigned>::const_iterator cellItr1 = hullNodeCells[*nodeItr1].begin();
             cellItr1 != hullNodeCells[*nodeItr1].end(); ++cellItr1){
           cellsOfInterest.insert(*cellItr1);

          // // Now the nodes of the hulls adjacent to our own. This second layer of adjacency
          // // catches neighbors on the full mesh that may not be neighbors on the hull mesh
          // // without adding a second communication step
          // for (set<unsigned>::const_iterator nodeItr2 = hullCellToNodes[*cellItr1].begin();
          //      nodeItr2 != hullCellToNodes[*cellItr1].end(); ++nodeItr2){
          //   for (set<unsigned>::const_iterator cellItr2 = hullNodeCells[*nodeItr2].begin();
          //        cellItr2 != hullNodeCells[*nodeItr2].end(); ++cellItr2){
          //     cellsOfInterest.insert(*cellItr2);
          //     // const unsigned otherProc = bisectSearch(domainCellOffset, *cellItr2);
          //     // if (otherProc != rank) neighborSet.insert(otherProc);
          //   }
          // }
        }
      }
    }

    // Now add the likely neighbors to our neighbor set
    for (set<unsigned>::const_iterator cellItr = cellsOfInterest.begin();
         cellItr != cellsOfInterest.end(); ++cellItr) {
      const unsigned otherProc = bisectSearch(domainCellOffset, *cellItr);
      if (otherProc != rank) neighborSet.insert(otherProc);
    }
    
    // // Now any hulls that have elements adjacent to ours in the hull mesh.
    // for (unsigned icell = domainCellOffset[rank]; icell != domainCellOffset[rank + 1]; ++icell) {
    //   for (typename vector<int>::const_iterator faceItr = hullMesh.cells[icell].begin();
    //        faceItr != hullMesh.cells[icell].end(); ++faceItr) {
    //     const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
    //     POLY_ASSERT(iface < hullMesh.faceCells.size());
    //     for (vector<unsigned>::const_iterator nodeItr = hullMesh.faces[iface].begin();
    //          nodeItr != hullMesh.faces[iface].end();
    //          ++nodeItr) {
    //       const unsigned inode = *nodeItr;
    //       POLY_ASSERT(inode < hullNodeCells.size());
    //       for (set<unsigned>::const_iterator itr = hullNodeCells[inode].begin();
    //            itr != hullNodeCells[inode].end();
    //            ++itr) {
    //         const unsigned otherProc = bisectSearch(domainCellOffset, *itr);
    //         if (otherProc != rank) neighborSet.insert(otherProc);
    //       }
    //     }
    //   }
    // }

    

    // Copy the neighbor set information the mesh.
    mesh.neighborDomains = vector<unsigned>();
    copy(neighborSet.begin(), neighborSet.end(), back_inserter(mesh.neighborDomains));
    sort(mesh.neighborDomains.begin(), mesh.neighborDomains.end());

    // // Blago!
    // cerr << "1st communicating with : ";
    // copy(mesh.neighborDomains.begin(), 
    //      mesh.neighborDomains.end(), 
    //      ostream_iterator<unsigned>(cerr, " "));
    // cerr << endl;
    // MPI_Barrier(MPI_COMM_WORLD);
    // // Blago!

#ifndef NDEBUG
    // Make sure everyone is consistent about who talks to whom.
    // This is potentially an expensive check on massively parallel systems!
    {
      for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
        unsigned numOthers = mesh.neighborDomains.size();
        vector<unsigned> otherNeighbors(mesh.neighborDomains);
        MPI_Bcast(&numOthers, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
        if (numOthers > 0) {
          otherNeighbors.resize(numOthers);
          MPI_Bcast(&(otherNeighbors.front()), numOthers, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
          POLY_ASSERT(rank == sendProc or
                 count(mesh.neighborDomains.begin(), mesh.neighborDomains.end(), sendProc) == 
                 count(otherNeighbors.begin(), otherNeighbors.end(), rank));
        }
      }
    }
#endif

    // For now we're not going to clever, just send all our generators to our 
    // potential neighbor set.  Pack 'em up.
    vector<char> localBuffer;
    serialize(generators, localBuffer);
    unsigned localBufferSize = localBuffer.size();

    // Fire off our sends asynchronously.
    vector<MPI_Request> sendRequests;
    sendRequests.reserve(2*mesh.neighborDomains.size());
    for (vector<unsigned>::const_iterator otherItr = mesh.neighborDomains.begin();
         otherItr != mesh.neighborDomains.end();
         ++otherItr) {
      const unsigned otherProc = *otherItr;
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&localBufferSize, 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &sendRequests.back());
      if (localBufferSize > 0) {
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&localBuffer.front(), localBufferSize, MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &sendRequests.back());
      }
    }
    POLY_ASSERT(sendRequests.size() <= 2*mesh.neighborDomains.size());

    // Get the info from each of our neighbors and append it to the result.
    // Simultaneously build the mapping of local generator ID to domain.
    for (vector<unsigned>::const_iterator otherItr = mesh.neighborDomains.begin();
         otherItr != mesh.neighborDomains.end();
         ++otherItr) {
      const unsigned otherProc = *otherItr;
      unsigned bufSize;
      MPI_Status recvStatus1, recvStatus2;
      MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 1, MPI_COMM_WORLD, &recvStatus1);
      if (bufSize > 0) {
        vector<char> buffer(bufSize, '\0');
        MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 2, MPI_COMM_WORLD, &recvStatus2);
        vector<char>::const_iterator itr = buffer.begin();
        vector<RealType> otherGenerators;
        deserialize(otherGenerators, itr, buffer.end());
        POLY_ASSERT(itr == buffer.end());
        copy(otherGenerators.begin(), otherGenerators.end(), back_inserter(generators));
        gen2domain.resize(generators.size()/Dimension, otherProc);
      }
    }

    // Make sure all our sends are completed.
    // cerr << rank << " : Waiting for generator sends to complete." << endl;
    vector<MPI_Status> sendStatus(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &sendRequests.front(), &sendStatus.front());
    // cerr << rank << " : DONE." << endl;
  }
  POLY_ASSERT(gen2domain.size() == generators.size()/Dimension);

  // // Denormalize the generator positions.
  // const unsigned ntotal = generators.size() / Dimension;
  // POLY_ASSERT(ntotal >= nlocal);
  // for (unsigned i = 0; i != ntotal; ++i) {
  //   for (unsigned j = 0; j != Dimension; ++j) {
  //     generators[Dimension*i + j] = generators[Dimension*i + j]/fscale + rlow[j];
  //   }
  // }
  

  // // Blago!
  // if (rank == 4) {
  //   cerr << "double points[" << generators.size() << "] = {" << endl;
  //   copy(generators.begin(), generators.end(), ostream_iterator<double>(cerr, ","));
  //   cerr << "};";
  // }
  // // Blago!


  // Construct the tessellation including the other domains' generators.
  this->tessellationWrapper(generators, mesh);


  // Blago!
#if HAVE_SILO
  if (siloBlago)
  {   
    POLY_ASSERT(gen2domain.size() == mesh.cells.size());
    vector<double> r2(mesh.cells.size());
    vector<double> px(mesh.cells.size());
    vector<double> py(mesh.cells.size());
    for (unsigned i = 0; i != mesh.cells.size(); ++i) {
      r2[i] = gen2domain[i];
      px[i] = generators[2*i];
      py[i] = generators[2*i+1];
    }
    map<string, double*> nodeFields, edgeFields, faceFields, cellFields;
    cellFields["domain"] = &r2[0];
    cellFields["gen_x"] = &px[0];
    cellFields["gen_y"] = &py[0];
    cerr << "Writing full mesh with " << mesh.cells.size() << endl;
    polytope::SiloWriter<Dimension, RealType>::write(mesh, nodeFields, edgeFields, faceFields, cellFields, "test_DistributedTessellator_fullMesh");
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif
  // Blago!

  // If requested, build the communication info for the shared nodes & faces
  // with our neighbor domains.
  if (mBuildCommunicationInfo and numProcs > 1) {

    // Build the reverse lookup from procID to index in the neighborDomain set.
    map<unsigned, unsigned> proc2offset;
    for (unsigned i = 0; i != mesh.neighborDomains.size(); ++i) proc2offset[mesh.neighborDomains[i]] = i;
    POLY_ASSERT(proc2offset.size() == mesh.neighborDomains.size());

    // Look for the faces we share with other processors.
    mesh.sharedFaces.resize(mesh.neighborDomains.size());
    for (int icell = 0; icell != nlocal; ++icell) {

      // Look for shared faces.
      const vector<int>& faces = mesh.cells[icell];
      for (vector<int>::const_iterator faceItr = faces.begin();
           faceItr != faces.end();
           ++faceItr) {
        const int iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
        POLY_ASSERT(mesh.faceCells[iface].size() == 1 or mesh.faceCells[iface].size() == 2);
        if (mesh.faceCells[iface].size() == 2) {
          const int jcell1 = mesh.faceCells[iface][0] < 0 ? ~mesh.faceCells[iface][0] : mesh.faceCells[iface][0];
          const int jcell2 = mesh.faceCells[iface][1] < 0 ? ~mesh.faceCells[iface][1] : mesh.faceCells[iface][1];
          POLY_ASSERT(jcell1 == icell or jcell2 == icell);
          const int jcell = jcell1 == icell ? jcell2 : jcell1;
          POLY_ASSERT(jcell < gen2domain.size());
          const unsigned otherProc = gen2domain[jcell];
          if (otherProc != rank) {
            POLY_ASSERT(proc2offset.find(otherProc) != proc2offset.end());
            const unsigned joff = proc2offset[otherProc];
            mesh.sharedFaces[joff].push_back(iface);
          }
        }
      }
    }

    // Look for shared nodes.
    mesh.sharedNodes.resize(mesh.neighborDomains.size());
    const vector<set<unsigned> > nodeCells = mesh.computeNodeCells();
    for (unsigned inode = 0; inode != nodeCells.size(); ++inode) {
      const set<unsigned>& cells = nodeCells[inode];
      for (typename set<unsigned>::const_iterator cellItr = cells.begin();
           cellItr != cells.end();
           ++cellItr) {
        POLY_ASSERT(*cellItr < gen2domain.size());
        const unsigned otherProc = gen2domain[*cellItr];
        if (otherProc != rank) {
          POLY_ASSERT(proc2offset.find(otherProc) != proc2offset.end());
          const unsigned joff = proc2offset[otherProc];
          mesh.sharedNodes[joff].push_back(inode);
        }
      }
    }

    
    // // Blago!
    // for (unsigned procID = 0; procID != numProcs; ++procID) {
    //   if (procID == rank) {
    //     for (unsigned iproc = 0; iproc != mesh.neighborDomains.size(); ++iproc) {
    //       cerr << rank << " <-> " << mesh.neighborDomains[iproc] << " : " 
    //            << mesh.nodes.size()/Dimension << endl;
    //       cerr << "  Nodes : ";
    //       copy(mesh.sharedNodes[iproc].begin(), mesh.sharedNodes[iproc].end(), ostream_iterator<unsigned>(cerr, " "));
    //       cerr << endl;
    //     }
    //   }
    //   MPI_Barrier(MPI_COMM_WORLD);
    // }
    // // Blago!



    // Remove any duplicate shared nodes.  (Faces should already be unique).
    for (unsigned i = 0; i != mesh.sharedNodes.size(); ++i) {
      sort(mesh.sharedNodes[i].begin(), mesh.sharedNodes[i].end());
      mesh.sharedNodes[i].erase(unique(mesh.sharedNodes[i].begin(), mesh.sharedNodes[i].end()), 
                                mesh.sharedNodes[i].end());
    }

    // // Blago!
    // for (unsigned iproc = 0; iproc != numProcs; ++iproc ){
    //    if (iproc==rank){
    //       cerr << endl << "PROC " << iproc << endl;
    //       for (unsigned i=0; i<mesh.nodes.size()/Dimension; ++i){
    //          cerr << "Node " << i << ":  ( " << mesh.nodes[Dimension*i] << " , " << mesh.nodes[Dimension*i+1] << " ) " << endl;
    //       }
    //    }
    //    MPI_Barrier(MPI_COMM_WORLD);
    // }
    // // Blago!


    // // Blago!
    // for (unsigned procID = 0; procID != numProcs; ++procID) {
    //   if (procID == rank) {
    //     for (unsigned iproc = 0; iproc != mesh.neighborDomains.size(); ++iproc) {
    //       cerr << rank << " <-> " << mesh.neighborDomains[iproc] << " : " 
    //            << mesh.nodes.size()/Dimension << endl;
    //       cerr << "  Nodes : ";
    //       copy(mesh.sharedNodes[iproc].begin(), mesh.sharedNodes[iproc].end(), ostream_iterator<unsigned>(cerr, " "));
    //       cerr << endl;
    //     }
    //   }
    //   MPI_Barrier(MPI_COMM_WORLD);
    // }
    // // Blago!

    // Compute the bounding box for the mesh coordinates.
    this->computeBoundingBox(mesh.nodes, rlow, rhigh);
    const RealType dx = DimensionTraits<Dimension, RealType>::maxLength(rlow, rhigh)/(1LL << 34);
    
    // Sort the shared elements.  This should make the ordering consistent on all domains without communication.
    unsigned numNeighbors = mesh.neighborDomains.size();
    for (unsigned idomain = 0; idomain != numNeighbors; ++idomain) {

      // Nodes.
      vector<Point> nodePoints;
      nodePoints.reserve(mesh.sharedNodes[idomain].size());
      for (vector<unsigned>::const_iterator itr = mesh.sharedNodes[idomain].begin();
           itr != mesh.sharedNodes[idomain].end();
           ++itr) nodePoints.push_back(DimensionTraits<Dimension, RealType>::constructPoint(&(mesh.nodes[Dimension * (*itr)]),
                                                                                            &rlow[0],
                                                                                            dx,
                                                                                            *itr));
      sort(nodePoints.begin(), nodePoints.end());
      for (unsigned i = 0; i != mesh.sharedNodes[idomain].size(); ++i) mesh.sharedNodes[idomain][i] = nodePoints[i].index;

      // Faces.
      vector<Point> facePoints;
      for (vector<unsigned>::const_iterator itr = mesh.sharedFaces[idomain].begin();
           itr != mesh.sharedFaces[idomain].end();
           ++itr) facePoints.push_back(DimensionTraits<Dimension, RealType>::faceCentroid(mesh,
                                                                                          *itr,
                                                                                          &rlow[0],
                                                                                          dx));
      sort(facePoints.begin(), facePoints.end());
      for (unsigned i = 0; i != mesh.sharedFaces[idomain].size(); ++i) mesh.sharedFaces[idomain][i] = facePoints[i].index;
    }
  }

  // // Blago!
  // for (unsigned procID = 0; procID != numProcs; ++procID) {
  //   if (procID == rank) {
  //     for (unsigned iproc = 0; iproc != mesh.neighborDomains.size(); ++iproc) {
  //       cerr << rank << " <-> " << mesh.neighborDomains[iproc] << " : " << mesh.nodes.size()/Dimension << endl;
  //       cerr << "  Nodes : ";
  //       copy(mesh.sharedNodes[iproc].begin(), mesh.sharedNodes[iproc].end(), ostream_iterator<unsigned>(cerr, " "));
  //       cerr << endl;
  //     }
  //   }
  //   MPI_Barrier(MPI_COMM_WORLD);
  // }
  // // Blago!

  // // Blago!
  // if (rank == 2) {
  //    vector<unsigned>::iterator it = 
  //       std::find(mesh.neighborDomains.begin(), mesh.neighborDomains.end(), 9);
  //    if (it != mesh.neighborDomains.end()) {
  //       unsigned badID = std::distance(mesh.neighborDomains.begin(), it);
  //       vector<Point2<double> > badNodes(mesh.sharedNodes[badID].size());
  //       int ii = 0;
  //       for (vector<unsigned>::const_iterator itr = mesh.sharedNodes[badID].begin();
  //            itr != mesh.sharedNodes[badID].end(); ++itr, ++ii) {
  //          badNodes[ii] = Point2<double>(mesh.nodes[2*(*itr)], mesh.nodes[2*(*itr)+1]);
  //       }
  //       cerr << "Shared nodes before deleteCells:\n";
  //       copy(mesh.sharedNodes[badID].begin(), mesh.sharedNodes[badID].end(), ostream_iterator<unsigned>(cerr, " "));
  //       cerr << endl;
  //       copy(badNodes.begin(), badNodes.end(), ostream_iterator<Point2<double> >(cerr, " "));
  //       cerr << endl;
  //    }
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
  // // Blago!


  // Remove the elements of the tessellation corresponding to
  // other domains generators, and renumber the resulting elements.
  vector<unsigned> cellMask(mesh.cells.size(), 1);
  fill(cellMask.begin() + nlocal, cellMask.end(), 0);
  deleteCells(mesh, cellMask);

  // Blago!
  // if (rank == 2) {
  //    vector<unsigned>::iterator it = 
  //       std::find(mesh.neighborDomains.begin(), mesh.neighborDomains.end(), 9);
  //    if (it != mesh.neighborDomains.end()) {
  //       unsigned badID = std::distance(mesh.neighborDomains.begin(), it);
  //       unsigned badNodeID = mesh.sharedNodes[badID][0];
  //       cerr << "Bad Node " << badNodeID << " @ (" << mesh.nodes[2*badNodeID]
  //            << "," << mesh.nodes[2*badNodeID+1] << ")" << endl;
  //    }
  // }
  // Blago!

  // // Blago!
  // for (unsigned iproc = 0; iproc != mesh.neighborDomains.size(); ++iproc) {
  //    cerr << rank << "\t<-->\t" << mesh.neighborDomains[iproc] 
  //         << "\t:\t" << mesh.sharedNodes[iproc].size()
  //         << "\t"    << mesh.sharedFaces[iproc].size() << endl;
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
  // // Blago!

  // Blago!
#if HAVE_SILO
  if (siloBlago)
  {
    vector<double> r2(mesh.cells.size());
    vector<double> px(mesh.cells.size());
    vector<double> py(mesh.cells.size());
    for (unsigned i = 0; i != mesh.cells.size(); ++i) {
       r2[i] = gen2domain[i];
       px[i] = generators[2*i];
       py[i] = generators[2*i+1];
    }
    map<string, double*> nodeFields, edgeFields, faceFields, cellFields;
    cellFields["domain"] = &r2[0];
    cellFields["gen_x"] = &px[0];
    cellFields["gen_y"] = &py[0];
    cerr << "Writing final mesh with " << mesh.cells.size() << endl;
    polytope::SiloWriter<Dimension, RealType>::write(mesh, nodeFields, edgeFields, faceFields, cellFields, "test_DistributedTessellator_finalMesh");
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif
  // Blago!


  // Remove any neighbors we don't actually share any info with.
  unsigned numNeighbors = mesh.neighborDomains.size();
  for (int i = numNeighbors - 1; i != -1; --i) {
    if (mesh.sharedNodes[i].size() == 0 and mesh.sharedFaces[i].size() == 0) {
      // cerr << "Removing neighbor " << i << " of " << numNeighbors << endl;
      mesh.neighborDomains.erase(mesh.neighborDomains.begin() + i);
      mesh.sharedNodes.erase(mesh.sharedNodes.begin() + i);
      mesh.sharedFaces.erase(mesh.sharedFaces.begin() + i);
    }
  }
  numNeighbors = mesh.neighborDomains.size();

  // // Blago!
  // cerr << "2nd communicating with : ";
  // copy(mesh.neighborDomains.begin(), 
  //      mesh.neighborDomains.end(), 
  //      ostream_iterator<unsigned>(cerr, " "));
  // cerr << endl;
  // MPI_Barrier(MPI_COMM_WORLD);
  // // Blago!



  // // Blago!
  // for (unsigned procID = 0; procID != numProcs; ++procID) {
  //   if (procID == rank) {
  //     for (unsigned iproc = 0; iproc != mesh.neighborDomains.size(); ++iproc) {
  //       cerr << rank << " <-> " << mesh.neighborDomains[iproc] << " : " << mesh.nodes.size()/Dimension << endl;
  //       cerr << "  Nodes : ";
  //       copy(mesh.sharedNodes[iproc].begin(), mesh.sharedNodes[iproc].end(), ostream_iterator<unsigned>(cerr, " "));
  //       cerr << endl;
  //     }
  //   }
  //   MPI_Barrier(MPI_COMM_WORLD);
  // }
  // // Blago!

  // In parallel we need to make sure the shared nodes are bit perfect the same.
  if (numProcs > 0) {

    // Figure out which domain owns the shared nodes.
    const unsigned nNodes = mesh.nodes.size()/Dimension;
    vector<unsigned> ownNode(nNodes, rank);
    for (unsigned idomain = 0; idomain != numNeighbors; ++idomain) {
      for (vector<unsigned>::const_iterator itr = mesh.sharedNodes[idomain].begin();
           itr != mesh.sharedNodes[idomain].end();
           ++itr) {
        POLY_ASSERT(*itr < nNodes);
        ownNode[*itr] = std::min(ownNode[*itr], mesh.neighborDomains[idomain]);
      }
    }

    // // Blago!
    // for (unsigned p = 0; p != numProcs; ++p) {
    //   if (p == rank) {
    //     for (unsigned i = 0; i != mesh.nodes.size()/Dimension; ++i) {
    //       if (abs(mesh.nodes[2*i] - 0.411765) + abs(mesh.nodes[2*i+1] - 0.235294) < 1.0e-5) {
    //         cerr << "Domain " << rank << " thinks node owned by " << ownNode[i] << endl;
    //       }
    //     }
    //   }
    //   MPI_Barrier(MPI_COMM_WORLD);
    // }
    // // Blago!

    // // Blago!
    // for (unsigned idomain = 0; idomain != numNeighbors; ++idomain){
    //    cerr << " shares " << mesh.sharedNodes[idomain].size() 
    //         << " nodes with " << mesh.neighborDomains[idomain] << endl;
    //    for (unsigned i=0; i != mesh.sharedNodes[idomain].size(); ++i){
    //       unsigned inode = mesh.sharedNodes[idomain][i];
    //       cerr << "   Node " << inode << " @ (" << mesh.nodes[2*inode] 
    //            << "," << mesh.nodes[2*inode+1] << ")" << endl;
    //    }
    // }      
    // // Blago!

    // Post the sends for any nodes we own, and note which nodes we expect to receive.
#ifndef NDEBUG
    vector<unsigned> bufSizes(numNeighbors, 0);
#endif
    list<vector<double> > sendCoords;
    list<vector<unsigned> > allRecvNodes;
    vector<MPI_Request> sendRequests;
    sendRequests.reserve(2*numNeighbors);
    for (unsigned idomain = 0; idomain != numNeighbors; ++idomain) {
      vector<unsigned> sendNodes;
      allRecvNodes.push_back(vector<unsigned>());
      vector<unsigned>& recvNodes = allRecvNodes.back();
      for (vector<unsigned>::const_iterator itr = mesh.sharedNodes[idomain].begin();
           itr != mesh.sharedNodes[idomain].end();
           ++itr) {
        POLY_ASSERT(ownNode[*itr] <= rank);
        if (ownNode[*itr] == rank) {
          sendNodes.push_back(*itr);
        } else if (ownNode[*itr] == mesh.neighborDomains[idomain]) {
          recvNodes.push_back(*itr);
        }
      }


      // // Blago!
      // for (unsigned ii=0; ii<sendNodes.size(); ++ii){
      //    unsigned inode = sendNodes[ii];
      //    cerr << "sending node " << inode << " @ (" << mesh.nodes[2*inode] << "," 
      //         << mesh.nodes[2*inode+1] << ") to processor " << mesh.neighborDomains[idomain] << endl;
      // }
      // for (unsigned ii=0; ii<recvNodes.size(); ++ii){
      //    unsigned inode = recvNodes[ii];
      //    cerr << "receiving node " << inode << " @ (" << mesh.nodes[2*inode] << "," 
      //         << mesh.nodes[2*inode+1] << ") from processor " << mesh.neighborDomains[idomain] << endl;
      // }
      // // Blago!


      // Send any nodes we have for this neighbor.
#ifndef NDEBUG
      bufSizes[idomain] = Dimension*sendNodes.size();
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&bufSizes[idomain], 1, MPI_UNSIGNED,
                mesh.neighborDomains[idomain], 9, MPI_COMM_WORLD, &sendRequests.back());
#endif
      if (sendNodes.size() > 0) {
        sendCoords.push_back(DimensionTraits<Dimension, RealType>::extractCoords(mesh.nodes, sendNodes));
        vector<RealType>& coords = sendCoords.back();
        POLY_ASSERT(coords.size() == Dimension*sendNodes.size());
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&coords.front(), Dimension*sendNodes.size(), DataTypeTraits<RealType>::MpiDataType(),
                  mesh.neighborDomains[idomain], 10, MPI_COMM_WORLD, &sendRequests.back());
      }
    }
    POLY_ASSERT(sendRequests.size() <= 2*numNeighbors);

    // Iterate over the neighbors again and look for any receive information.
    list<vector<unsigned> >::const_iterator recvNodesItr = allRecvNodes.begin();
    for (unsigned idomain = 0; idomain != numNeighbors; ++idomain, ++recvNodesItr) {
      POLY_ASSERT(recvNodesItr != allRecvNodes.end());
      const vector<unsigned>& recvNodes = *recvNodesItr;

#ifndef NDEBUG      
      unsigned otherSize;
      MPI_Status recvStatus;
      MPI_Recv(&otherSize, 1, MPI_UNSIGNED, mesh.neighborDomains[idomain], 9, MPI_COMM_WORLD, &recvStatus);
      POLY_ASSERT2(otherSize == Dimension*recvNodes.size(),
              "Bad message size (" << mesh.neighborDomains[idomain] << " " << otherSize << ") (" << rank << " " << Dimension*recvNodes.size() << ")");
#endif

      if (recvNodes.size() > 0) {
        vector<RealType> recvCoords(Dimension*recvNodes.size());
        MPI_Status recvStatus;
        MPI_Recv(&recvCoords.front(), Dimension*recvNodes.size(), DataTypeTraits<RealType>::MpiDataType(),
                 mesh.neighborDomains[idomain], 10, MPI_COMM_WORLD, &recvStatus);

        // // Blago!
        // if (otherSize != Dimension*recvNodes.size()) {
        //   cout << "Bad message size from " << mesh.neighborDomains[idomain] << " being sent to " << rank << endl;
        //   cout << "Other coordinates : ";
        //   copy(recvCoords.begin(), recvCoords.end(), ostream_iterator<RealType>(cout, " "));
        //   cout << endl;
        //   cout << "My coordinates : ";
        //   for (unsigned j = 0; j != recvNodes.size(); ++j) {
        //     const unsigned i = recvNodes[j];
        //     for (unsigned k = 0; k != Dimension; ++k) cout << mesh.nodes[Dimension*i + k] << " ";
        //   }
        //   cout << endl;
        //   POLY_ASSERT(false);
        // }
        // // Blago!

        // Unpack the coordinates to the receive nodes.
        for (unsigned j = 0; j != recvNodes.size(); ++j) {
          const unsigned i = recvNodes[j];
          for (unsigned k = 0; k != Dimension; ++k) mesh.nodes[Dimension*i + k] = recvCoords[Dimension*j + k];
        }
      }
    }

    // Wait until all our send are complete.
    vector<MPI_Status> sendStatus(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &sendRequests.front(), &sendStatus.front());
  }

  // Post-conditions.
#ifndef NDEBUG
  const string msg = checkDistributedTessellation(mesh);
  if (msg != "ok" and rank == 0) cerr << msg;
  //POLY_ASSERT(msg == "ok");
#endif
}

//------------------------------------------------------------------------------
// Dispatch a serial tessellation call appropriately.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
tessellationWrapper(const vector<RealType>& points,
                    Tessellation<Dimension, RealType>& mesh) const {
  switch (mType) {
  case unbounded:
    mSerialTessellator->tessellate(points, mesh);
    break;

  case box:
    POLY_ASSERT(mLow != 0);
    POLY_ASSERT(mHigh != 0);
    mSerialTessellator->tessellate(points, mLow, mHigh, mesh);
    break;

  case plc:
    POLY_ASSERT(mPLCptr != 0);
    mSerialTessellator->tessellate(points, *mPLCpointsPtr, *mPLCptr, mesh);
    break;
  }
}

//------------------------------------------------------------------------------
// Compute the bouding box to encompass all the points.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
DistributedTessellator<Dimension, RealType>::
computeBoundingBox(const vector<RealType>& points,
                   RealType* rlow,
                   RealType* rhigh) const {

  for (unsigned i = 0; i != Dimension; ++i) {
    rlow[i] = numeric_limits<RealType>::max();
    rhigh[i] = (numeric_limits<RealType>::is_signed ? -rlow[i] : numeric_limits<RealType>::min());
  }

  // Find the local min & max.
  switch (mType) {
  case unbounded:
    geometry::computeBoundingBox<Dimension, RealType>(points, false, rlow, rhigh);
    break;

  case box:
    POLY_ASSERT(mLow != 0);
    POLY_ASSERT(mHigh != 0);
    for (unsigned j = 0; j != Dimension; ++j) {
      rlow[j] = mLow[j];
      rhigh[j] = mHigh[j];
    }
    return;

  case plc:
    POLY_ASSERT(mPLCptr != 0);
    POLY_ASSERT(mPLCpointsPtr != 0);
    for (vector<vector<int> >::const_iterator facetItr = mPLCptr->facets.begin();
         facetItr != mPLCptr->facets.end();
         ++facetItr) {
      for (vector<int>::const_iterator iItr = facetItr->begin();
           iItr != facetItr->end();
           ++iItr) {
        const unsigned i = *iItr;
        for (unsigned j = 0; j != Dimension; ++j) {
          // rlow[j] =  min(rlow[j],  points[Dimension*i + j]);
          // rhigh[j] = max(rhigh[j], points[Dimension*i + j]);
          rlow[j]  = min(rlow[j],  (*mPLCpointsPtr)[Dimension*i + j]);
          rhigh[j] = max(rhigh[j], (*mPLCpointsPtr)[Dimension*i + j]);
        }
      }
    }

    break;
  }

  // Find the global results.
  for (unsigned j = 0; j != Dimension; ++j) {
    rlow[j] = allReduce(rlow[j], MPI_MIN, MPI_COMM_WORLD);
    rhigh[j] = allReduce(rhigh[j], MPI_MAX, MPI_COMM_WORLD);
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class DistributedTessellator<2, double>;
template class DistributedTessellator<3, double>;

}
