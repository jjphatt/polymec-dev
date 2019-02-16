// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BLOCKMESH_H
#define POLYMEC_BLOCKMESH_H

#include "geometry/coord_mapping.h"
#include "geometry/unimesh.h"

/// \addtogroup geometry geometry
///@{

/// \class blockmesh
/// A blockmesh, or block mesh, is a collection of three-dimensional cartesian 
/// meshes (unimeshes) called blocks. These blocks can be connected to one 
/// another as long as their boundaries have compatible dimensions (numbers of
/// cells). The block mesh manages these blocks and their connectivity.
///
/// A blockmesh is distributed across processes in an MPI communicator. Every 
/// process contains all of the blocks in a blockmesh. Each of these blocks
/// is domain-decomposed across these processes.
///
/// Blocks in a block mesh are stitched together by identifying nodes on the 
/// common face between two blocks. Each block has 8 nodes. 
///
/// Looking down (in the -z direction) on a block in the xy plane, nodes 0-3 
/// traverse the bottom block face counterclockwise, starting with the "lower 
/// left" node. Nodes 4-7 traverse the top block face in the same way. This 
/// is the standard way that nodes are indexed in linear hexahedral elements 
/// in the finite element method.
///
/// There are no periodic blocks in a block mesh. If you want a periodic block, 
/// you can always connectg the opposite faces of that block to get the same
/// effect.
typedef struct blockmesh_t blockmesh_t;

//------------------------------------------------------------------------
//                          Construction methods
//------------------------------------------------------------------------
// The following methods are used to construct block meshs.
// blockmesh_finalize() must be called after a mesh has been properly
// constructed.
//------------------------------------------------------------------------

/// Creates a new blockmesh with no blocks, defined on the given communicator.
/// The blockmesh can house unimesh blocks whose patches have the given dimensions.
/// \param [in] comm The MPI communicator for the blockmesh.
/// \param [in] patch_nx The number of "x" cells in a patch within each block of the blockmesh.
/// \param [in] patch_ny The number of "y" cells in a patch within each block of the blockmesh.
/// \param [in] patch_nx The number of "z" cells in a patch within each block of the blockmesh.
/// \memberof blockmesh
blockmesh_t* blockmesh_new(MPI_Comm comm, int patch_nx, int patch_ny, int patch_nz);

/// Adds a new empty block to the blockmesh with the given numbers of patches in the "x", "y", 
/// and "z" directions. The block is empty in the sense that it contains no patches.
/// \param [in] block_domain A bounding box containing the "logical" domain of the block's
///                          coordinate mapping (e.g. [-1, +1] x [-1, +1] x [-1, +1]).
/// \param [in] block_coords A coordinate mapping from the block's domain to the new block's 
///                          local coordinate system. 
///                          * Must be non-NULL, since the coordinate systems for the blocks 
///                            within the mesh must form a smooth atlas of compatible charts 
///                            (permitting diffeomorphisms between charts).
///                          * This diffeomorphism also requires block_coords to have an inverse.
///                          * The blockmesh assumes ownership of block_coords, so you must 
///                            retain a reference to continue using it outside of this context.
/// \param [in] num_x_patches The number of patches in the "x" direction within the new block.
/// \param [in] num_y_patches The number of patches in the "y" direction within the new block.
/// \param [in] num_z_patches The number of patches in the "z" direction within the new block.
/// \returns the index of the new block within the blockmesh.
/// \memberof blockmesh
int blockmesh_add_block(blockmesh_t* mesh, 
                        bbox_t* block_domain,
                        coord_mapping_t* block_coords,
                        int num_x_patches, 
                        int num_y_patches, 
                        int num_z_patches);

/// Returns the block boundary associated with the given set of block nodes, or -1 
/// if the nodes do not match any of the block's faces. 
/// \param [in] block_nodes An array of 4 block nodes that supposedly match those
///                         in one of the block's 6 faces. The order in which 
///                         the nodes are specified doesn't matter.
/// \returns An integer between 0 and 5 that can be cast to a unimesh_boundary_t, or -1 
///          if the given nodes don't match with those on a block boundary.
/// \memberof blockmesh
int blockmesh_block_boundary_for_nodes(blockmesh_t* mesh, int block_nodes[4]);

/// Returns true if the two blocks with the given indices within the blockmesh
/// can be connected successfully, false if not.
/// \param [in] block1_index The index of the first of the two blocks to be 
///                          connected within the mesh.
/// \param [in] block1_nodes An array containing the 4 nodes in the first block
///                          to be identified with the corresponding nodes 
///                          in the second block (block2_nodes).
/// \param [in] block2_index The index of the second of the two blocks to be 
///                          connected within the mesh.
/// \param [in] block2_nodes An array containing the 4 nodes in the second block
///                          to be identified with the corresponding nodes 
///                          in the first block (block1_nodes).
/// \param [out] reason If non-NULL, this string stores an internal string 
///                     that describes any condition preventing a successful
///                     block connection.
/// \returns True if the two blocks can be connected, false if not.
/// \memberof blockmesh
bool blockmesh_can_connect_blocks(blockmesh_t* mesh, 
                                  int block1_index, 
                                  int block1_nodes[4],
                                  int block2_index,
                                  int block2_nodes[4],
                                  char** reason);

/// Connects two blocks with the given indices within a block mesh in a manner 
/// specified by parameters. You must call this function on every process 
/// within the mesh's communicator, and the blocks must be connected in such a 
/// way that all block dimensions are compatible. The two indices can refer to 
/// the same block, in which case the block connects to itself.
/// \note This function assumes that the blocks can be successfully connected.
/// \param [in] block1_index The index of the first of the two blocks to be 
///                          connected within the mesh.
/// \param [in] block1_nodes An array containing the 4 nodes in the first block
///                          to be identified with the corresponding nodes 
///                          in the second block (block2_nodes).
/// \param [in] block2_index The index of the second of the two blocks to be 
///                          connected within the mesh.
/// \param [in] block2_nodes An array containing the 4 nodes in the second block
///                          to be identified with the corresponding nodes 
///                          in the first block (block1_nodes).
/// \memberof blockmesh
/// \collective
void blockmesh_connect_blocks(blockmesh_t* mesh, 
                              int block1_index, 
                              int block1_nodes[4],
                              int block2_index,
                              int block2_nodes[4]);

/// Adds patches to all blocks within the block mesh, assigning them to 
/// processes linearly according to a naive partitioning. This is really only
/// intended to get you off the ground, and not to provide a load-balanced 
/// partitioning. After calling this to assign the patches, call 
/// \ref blockmesh_finalize and then use \ref repartition_blockmesh to do load 
/// balancing.
void blockmesh_assign_patches(blockmesh_t* mesh);

/// Finalizes the construction process for the block mesh. The function checks 
/// the following conditions, throwing a fatal error if any are not met:
/// * Every patch in every block in the mesh must have been inserted on some 
///   process in the communicator.
/// * Every block-to-block connection must be consistent.
/// This function must be called before any of the mesh's usage methods are 
/// invoked. It should only be called once.
/// \memberof blockmesh
void blockmesh_finalize(blockmesh_t* mesh);

//------------------------------------------------------------------------
//                          Usage methods
//------------------------------------------------------------------------
// The following methods can only be used after a blockmesh has been 
// fully constructed and finalized.
//------------------------------------------------------------------------

/// Returns true if the given mesh has been finalized, false otherwise.
/// \memberof blockmesh
bool blockmesh_is_finalized(blockmesh_t* mesh);

/// Destroys the given mesh and all of its patches.
/// \memberof blockmesh
void blockmesh_free(blockmesh_t* mesh);

/// Returns the MPI communicator on which the blockmesh is defined.
/// \memberof blockmesh
MPI_Comm blockmesh_comm(blockmesh_t* mesh);

/// Returns the number of blocks within the mesh.
/// \memberof blockmesh
int blockmesh_num_blocks(blockmesh_t* mesh);

/// Returns the block with the given index within the mesh.
/// \param [in] index The index of the requested block.
/// \memberof blockmesh
unimesh_t* blockmesh_block(blockmesh_t* mesh, int index);

/// Returns true if the block with the given index in the mesh is connected
/// to another block in the mesh on the given boundary, false if not.
/// \param [in] index The index of the block.
/// \param [in] boundary The boundary of the block for the connection in question.
/// \memberof blockmesh
bool blockmesh_block_is_connected(blockmesh_t* field,
                                  int index,
                                  unimesh_boundary_t boundary);

/// Allows the traversal of all blocks in the blockmesh. 
/// \param [inout] pos Stores the index of the next block in the mesh. 
///                    Set *pos to 0 to reset the traversal.
/// \param [out] block Stores the next block in the mesh.
/// \param [out] block_domain If not NULL, stores the bounding box 
///                           representing the logical domain of the 
///                           next block.
/// \param [out] block_coords If not NULL, stores the coordinate 
///                           mapping for the representing the next 
///                           block.
/// \returns True if the mesh contains another block, false if not.
bool blockmesh_next_block(blockmesh_t* mesh, 
                          int* pos, 
                          unimesh_t** block,
                          bbox_t* block_domain,
                          coord_mapping_t** block_coords);

typedef struct blockmesh_field_t blockmesh_field_t;

/// Repartitions the given blockmesh and redistributes data to each of the 
/// given fields. 
/// \param [inout] mesh A pointer that stores the old mesh, which is consumed and 
///                     replaced with the repartitioned mesh.
/// \param [in] weights If non-NULL, this is an array containing an integer weight
///                     for each locally-stored patch within the blockmesh. The 
///                     weights can be assigned with a nested traversal of 
///                     blocks over the mesh, and patches over each block.
/// \param [in] imbalance_tol A tolerance that governs the partitioning. The load
///                           imbalance produced by the repartitioning doesn't 
///                           exceed this tolerance.
/// \param [inout] fields An array of field pointers containing fields whose data
///                       is repartitioned in tandem with the blockmesh.
/// \param [in] num_fields The length of the fields array.
/// \note In addition, each repartitioned field needs to have any boundary 
/// conditions reinstated, since these boundary conditions are not 
/// transmitted between processes.
/// \relates blockmesh
/// \collective Collective on mesh's communicator.
void repartition_blockmesh(blockmesh_t** mesh, 
                           int* weights,
                           real_t imbalance_tol,
                           blockmesh_field_t** fields,
                           size_t num_fields);

///@}

#endif

