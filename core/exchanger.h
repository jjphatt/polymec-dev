// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_EXCHANGER_H
#define POLYMEC_EXCHANGER_H

#include "core/polymec.h"
#include "core/unordered_map.h"
#include "core/serializer.h"

/// \addtogroup core core
///@{

/// \class exchanger
/// This type implements an MPI transmitter/receiver for exchanging 
/// data between processes in a point-to-point fashion. Objects of this type 
/// are garbage collected.
typedef struct exchanger_t exchanger_t;

/// Constructs a new exchanger on the given communicator.
/// \memberof exchanger
exchanger_t* exchanger_new(MPI_Comm comm);

/// Constructs a new exchanger on the given communicator with the given rank.
/// \memberof exchanger
exchanger_t* exchanger_new_with_rank(MPI_Comm comm, int rank);

/// Creates a complete copy of the given exchanger.
/// \memberof exchanger
exchanger_t* exchanger_clone(exchanger_t* ex);

/// Returns the MPI communicator on which this exchanger is defined.
/// \memberof exchanger
MPI_Comm exchanger_comm(exchanger_t* ex);

/// Exchanges data of the given type in the given array with other processors.
/// \param [in,out] data An array of data for which values is exchanged with (sent and received to
///                      and from) other processes by this exchanger.
/// \param [in] stride The stride (number of elements) of data exchanged at each index of the array.
/// \param [in] tag The MPI tag used in the underlying exchange of data.
/// \param [in] type The MPI type for the data being exchanged.
/// \memberof exchanger
void exchanger_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type);

/// Begins an asynchronous data exchange. 
/// \param [in] data An array of data for which values is exchanged with (sent and received to
///                  and from) other processes by this exchanger.
/// \param [in] stride The stride (number of elements) of data exchanged at each index of the array.
/// \param [in] tag The MPI tag used in the underlying exchange of data.
/// \param [in] type The MPI type for the data being exchanged.
/// \returns a unique token that can be passed to \ref exchanger_finish_exchange to finish the exchange.
/// \memberof exchanger
int exchanger_start_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type);

/// Concludes the asynchronous exchange corresponding to the given token.
/// This fills the array given in \ref exchanger_start_exchange with the data it expects.
/// \param [in] token A token returned by \ref exchanger_start_exchange.
/// \memberof exchanger
void exchanger_finish_exchange(exchanger_t* ex, int token);

/// Returns the maximum index to be sent by this exchanger.
/// \memberof exchanger
int exchanger_max_send(exchanger_t* ex);

/// Returns the maximum index to be received by this exchanger.
/// \memberof exchanger
int exchanger_max_receive(exchanger_t* ex);

/// Enables deadlock detection, setting the threshold to the given number of 
/// seconds. Deadlocks will be reported to the given rank on the given stream.
/// \param [in] threshold The number of seconds after which an exchanger transaction 
///                       (a send or receive) is considered to have hung because of a 
///                       deadlock condition.
/// \param [in] output_rank The MPI rank on which any diagnostic output is reported
///                         for a deadlock condition.
/// \param [in] stream If non-NULL on the given output rank, specifies the stream to which 
///                    diagnostic output is written.
/// \memberof exchanger
void exchanger_enable_deadlock_detection(exchanger_t* ex, 
                                         real_t threshold,
                                         int output_rank,
                                         FILE* stream);

/// Disables deadlock detection.
/// \memberof exchanger
void exchanger_disable_deadlock_detection(exchanger_t* ex);

/// Returns true if deadlock detection is enabled, false otherwise.
/// \memberof exchanger
bool exchanger_deadlock_detection_enabled(exchanger_t* ex);

/// This writes a string representation of the exchanger to the given file stream.
/// \memberof exchanger
void exchanger_fprintf(exchanger_t* ex, FILE* stream);

/// This creates a serializer object that can read and write exchangers to byte streams.
/// \memberof exchanger
serializer_t* exchanger_serializer(void);

/// Establishes a communication pattern in which this exchanger sends data at 
/// the given indices of an array to the given remote process. Note that 
/// remote_process must differ from the local rank on the exchanger's communicator.
/// \param [in] remote_process The rank of the process to which the data at the given indices is 
///                            sent. Must differ from the local rank on the exchanger's communicator.
/// \param [in] indices The indices identifying locations in a data array to be sent to the remote process.
/// \param [in] num_indices The number of indices in the indices array.
/// \param [in] copy_indices If true, indices will be copied to this exchanger. If false, the exchanger 
///                          assumes ownership of the indices array passed to this function.
/// \memberof exchanger
void exchanger_set_send(exchanger_t* ex, 
                        int remote_process, 
                        int* indices, 
                        int num_indices, 
                        bool copy_indices);

/// Establishes communications patterns in which this exchanger sends data at 
/// the given indices of an array to various remote processes. 
/// \param [in] send_map A mapping of remote process ranks to \ref int_array objects 
///                      containing local indices identifying data to be sent. 
///                      Each remote_process must differ from the local rank on the 
///                      exchanger's communicator.
/// \memberof exchanger
void exchanger_set_sends(exchanger_t* ex, int_ptr_unordered_map_t* send_map);

/// Sets an offset for indices in data arrays that are sent to other processes.
/// This can be used to allow multiple exchangers to exchange data correctly 
/// in arrays that are aggregates of data associated with different distributed 
/// objects.
/// \memberof exchanger
void exchanger_set_send_offset(exchanger_t* ex, ssize_t offset);

/// Returns the number of processes to which this exchanger sends data.
/// \memberof exchanger
int exchanger_num_sends(exchanger_t* ex);

/// Removes the given remote process from the set of processes to which this 
/// exchanger sends data.
/// \memberof exchanger
void exchanger_delete_send(exchanger_t* ex, int remote_process);

/// Allows the traversal of the set of send indices for remote processes.
/// \param pos [in,out] Controls the traversal. Set to 0 to reset.
/// \param remote_process [out] Stores the remote process in the next send transaction.
/// \param indices [out] Stores an internal pointer to the array of indices identifying data 
///                      sent to the remote process.
/// \param num_indices [out] Stores the length of the indices array.
/// \returns true if another send transaction is available in the exchanger, false otherwise.
/// \memberof exchanger
bool exchanger_next_send(exchanger_t* ex, 
                         int* pos, 
                         int* remote_process, 
                         int** indices, 
                         int* num_indices);

/// Retrieves the indices and number of indices for a send transaction 
/// with the given process.
/// \param remote_process [in] The remote process in the specified send transaction.
/// \param indices [out] Stores an internal pointer to the array of indices identifying data 
///                      sent to the remote process.
/// \param num_indices [out] Stores the length of the indices array.
/// \returns true if such a transaction exists, false if not.
/// \memberof exchanger
bool exchanger_get_send(exchanger_t* ex, int remote_process, int** indices, int* num_indices);

/// Establishes a communication pattern in which this exchanger receives data at 
/// the given indices of an array from the given remote process. Note that
/// \param [in] remote_process The rank of the process from which the data at the given indices is 
///                            received. Must differ from the local rank on the exchanger's communicator.
/// \param [in] indices The indices identifying locations in a data array to be received from the remote 
///                     process.
/// \param [in] num_indices The number of indices in the indices array.
/// \param [in] copy_indices If true, indices will be copied to this exchanger. If false, the exchanger 
///                          assumes ownership of the indices array passed to this function.
/// \memberof exchanger
void exchanger_set_receive(exchanger_t* ex, 
                           int remote_process, 
                           int* indices, 
                           int num_indices, 
                           bool copy_indices);

/// Establishes communications patterns in which this exchanger receives data at 
/// the given indices of an array from various remote processes. Here, recv_map 
/// maps remote process ranks to int_arrays containing local indices identifying
/// locations where received data will be stored. Note that the remote_processes must differ from the 
/// local rank on the exchanger's communicator.
/// \param [in] recv_map A mapping of remote process ranks to \ref int_array objects 
///                      containing local indices identifying data to be received. 
///                      Each remote_process must differ from the local rank on the 
///                      exchanger's communicator.
/// \memberof exchanger
void exchanger_set_receives(exchanger_t* ex, int_ptr_unordered_map_t* recv_map);

/// Sets an offset for indices in data arrays that are received from other 
/// processes. This can be used to allow multiple exchangers to exchange data 
/// correctly in arrays that are aggregates of data associated with different 
/// distributed objects.
/// \memberof exchanger
void exchanger_set_receive_offset(exchanger_t* ex, ssize_t offset);

/// Returns the number of processes from which this exchanger receives data.
/// \memberof exchanger
int exchanger_num_receives(exchanger_t* ex);

/// Removes the given remote process from the set of processes from which this 
/// exchanger receives data.
void exchanger_delete_receive(exchanger_t* ex, int remote_process);

/// Allows the traversal of the set of receive indices for remote processes.
/// \param pos [in,out] Controls the traversal. Set to 0 to reset.
/// \param remote_process [out] Stores the remote process in the next receive transaction.
/// \param indices [out] Stores an internal pointer to the array of indices identifying data 
///                      received from the remote process.
/// \param num_indices [out] Stores the length of the indices array.
/// \returns true if another receive transaction is available in the exchanger, false otherwise.
/// \memberof exchanger
bool exchanger_next_receive(exchanger_t* ex, 
                            int* pos, 
                            int* remote_process, 
                            int** indices, 
                            int* num_indices);

/// Retrieves the indices and number of indices for a receive transaction 
/// with the given process. 
/// \param remote_process [in] The remote process in the specified receive transaction.
/// \param indices [out] Stores an internal pointer to the array of indices identifying data 
///                      sent to the remote process.
/// \param num_indices [out] Stores the length of the indices array.
/// \returns true if such a transaction exists, false if not.
/// \memberof exchanger
bool exchanger_get_receive(exchanger_t* ex, int remote_process, int** indices, int* num_indices);

/// Verifies the consistency of the exchanger.
/// This function is expensive and involves parallel communication. It must be called 
/// by all processes on the communicator for the exchanger. 
/// \param [in] handler A function that accepts a formatted string. If non-NULL, 
/// this function is called with a string describing any errors encountered. 
/// \returns true if the verification succeeds, false if not. 
/// \memberof exchanger
/// \collective Collective on the exchanger's communicator.
bool exchanger_verify(exchanger_t* ex, void (*handler)(const char* format, ...));

//------------------------------------------------------------------------
//                      Exchanging parallel metadata
//------------------------------------------------------------------------
// It is often expedient to send metadata between communicating processes
// instead of performing an exchange on full data arrays. The following 
// methods allow the transfer of metadata from a set of "send arrays" 
// to a set of "receive arrays", based on the topology established within 
// the given exchanger.
//------------------------------------------------------------------------

/// Creates a set of arrays to use in sending metadata of the given type. This 
/// allocates and returns a set of exchanger_num_sends(ex) arrays, each sized 
/// and indexed appropriately for sending metadata with the given stride using 
/// the given exchanger.
/// \memberof exchanger
void** exchanger_create_metadata_send_arrays(exchanger_t* ex,
                                             MPI_Datatype type,
                                             int stride);

/// Frees the storage associated with the metadata send arrays created with 
/// exchanger_create_metadata_send_arrays.
/// \memberof exchanger
void exchanger_free_metadata_send_arrays(exchanger_t* ex,
                                         void** arrays);

/// Create a set of arrays to use for receiving metadata of the given type. 
/// This allocates and returns a set of exchanger_num_receives(ex) arrays, each 
/// sized and indexed appropriately for receiving metadata with the given stride 
/// using the given exchanger.
/// \memberof exchanger
void** exchanger_create_metadata_receive_arrays(exchanger_t* ex,
                                                MPI_Datatype type,
                                                int stride);

/// Frees the storage associated with the metadata receive arrays created with 
/// exchanger_create_metadata_receive_arrays.
/// \memberof exchanger
void exchanger_free_metadata_receive_arrays(exchanger_t* ex,
                                            void** arrays);
 
/// \enum exchanger_metadata_dir_t
/// This enumerated type indicates the direction of metadata transfer. 
/// EX_METADATA_FORWARD indicates a "forward" transfer (from sends to receives), 
/// while EX_METADATA_REVERSE indicates a "reverse" transer (from receives to sends).
typedef enum
{
  EX_METADATA_FORWARD,
  EX_METADATA_REVERSE
} exchanger_metadata_dir_t;

/// Performs a transfer of metadata of the given type, with the given stride, 
/// in the given direction.
/// \memberof exchanger
void exchanger_transfer_metadata(exchanger_t* ex,
                                 void** send_arrays,
                                 void** receive_arrays,
                                 int stride,
                                 int tag,
                                 MPI_Datatype type,
                                 exchanger_metadata_dir_t direction);

/// Starts an asyncronous transfer of metadata, returning a token that can be used 
/// to finish it.
/// \memberof exchanger
int exchanger_start_metadata_transfer(exchanger_t* ex,
                                      void** send_arrays,
                                      void** receive_arrays,
                                      int stride,
                                      int tag,
                                      MPI_Datatype type,
                                      exchanger_metadata_dir_t direction);

/// Finishes the asyncronous metadata transfer associated with the given token.
/// \memberof exchanger
void exchanger_finish_metadata_transfer(exchanger_t* ex,
                                        int token);

///@}

#endif
