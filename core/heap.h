// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_HEAP_H
#define POLYMEC_HEAP_H

#include "core/comparators.h"
#include "core/polymec.h"

// A heap is a data structure that keeps its contents sorted.
// One defines a heap using
// DEFINE_HEAP(heap_name, element, comparator)
//
// Interface for a type x_heap_t (with element type x) defined with 
// DEFINE_HEAP(x_heap, x, x_comparator):
// 
// x_heap_t* x_heap_new() - Creates a new empty heap.
// x_heap_t* x_heap_with_comparator(comp) - Creates a new empty heap that uses the given comparator.
// x_heap_t* x_heap_with_reversed_comparator(comp) - Creates a new empty heap that uses the given comparator, reversed.
// void x_heap_free(x_heap_t* heap) - Destroys the heap.
// void x_heap_push_with_dtor(x_heap_t* heap, x datum) - Inserts a datum into the heap with a destructor.
// void x_heap_push(x_heap_t* heap, x datum) - Inserts a datum into the heap.
// void x_heap_pop(x_heap_t* heap) - Removes the largest element from the heap.
// x x_heap_front(x_heap_t* heap) - Returns the largest element in the heap.
// bool x_heap_empty(x_heap_t* heap) - Returns true if the heap is empty, false otherwise.
// void x_heap_clear(x_heap_t* heap) - Clears the contents of the heap.
#define DEFINE_HEAP(heap_name, element, comparator) \
typedef element heap_name##_element_t; \
typedef struct heap_name##_t heap_name##_t; \
typedef int (*heap_name##_comparator_t)(element left, element right); \
typedef void (*heap_name##_dtor)(element e); \
struct heap_name##_t \
{ \
	int capacity; \
	int size; \
	element* data; \
  heap_name##_dtor* dtors; \
  heap_name##_comparator_t cmp; \
  int cmp_factor; \
}; \
\
static inline heap_name##_t* heap_name##_with_comparator(heap_name##_comparator_t comp) \
{ \
  heap_name##_t* heap = polymec_malloc(sizeof(heap_name##_t)); \
  heap->capacity = 4; \
  heap->size = 0; \
  heap->data = polymec_malloc(sizeof(element) * heap->capacity); \
  heap->cmp = comp; \
  heap->cmp_factor = 1; \
  return heap; \
} \
\
static inline heap_name##_t* heap_name##_new() \
{ \
  return heap_name##_with_comparator(comparator); \
} \
\
static inline heap_name##_t* heap_name##_with_reversed_comparator(heap_name##_comparator_t comp) \
{ \
  heap_name##_t* heap = heap_name##_with_comparator(comp); \
  heap->cmp_factor = -1; \
  return heap; \
} \
\
static inline void heap_name##_free(heap_name##_t* heap) \
{ \
  for (int i = 0; i < heap->size; ++i) \
  { \
    if (heap->dtors[i] != NULL) \
    { \
      heap->dtors[i](heap->data[i]); \
    } \
  } \
  polymec_free(heap->data); \
  polymec_free(heap->dtors); \
  polymec_free(heap); \
} \
\
static inline void heap_name##_push_with_dtor(heap_name##_t* heap, heap_name##_element_t value, heap_name##_dtor dtor) \
{ \
	if (heap->size == heap->capacity) \
	{ \
		heap->capacity *= 2; \
		heap->data = polymec_realloc(heap->data, sizeof(element) * heap->capacity); \
		heap->dtors = polymec_realloc(heap->dtors, sizeof(heap_name##_dtor) * heap->capacity); \
	} \
  int index, parent; \
	for(index = heap->size++; index; index = parent) \
	{ \
		parent = (index - 1) >> 1; \
		if (heap->cmp_factor * heap->cmp(heap->data[parent], value) > 0) break; \
		heap->data[index] = heap->data[parent]; \
		heap->dtors[index] = heap->dtors[parent]; \
	} \
	heap->data[index] = value; \
	heap->dtors[index] = dtor; \
} \
\
static inline void heap_name##_push(heap_name##_t* heap, heap_name##_element_t value) \
{ \
  heap_name##_push_with_dtor(heap, value, NULL); \
} \
\
static inline void heap_name##_pop(heap_name##_t* heap) \
{ \
  --heap->size; \
	element temp = heap->data[heap->size]; \
	heap_name##_dtor temp_dtor = heap->dtors[heap->size]; \
	if ((heap->size <= (heap->capacity >> 2)) && (heap->capacity > 4)) \
	{ \
		heap->capacity >>= 1; \
		heap->data = polymec_realloc(heap->data, sizeof(element) * heap->capacity); \
		heap->dtors = polymec_realloc(heap->dtors, sizeof(heap_name##_dtor) * heap->capacity); \
	} \
	int index, swap; \
	for(index = 0; 1; index = swap) \
	{ \
		swap = (index << 1) + 1; \
		if (swap >= heap->size) break; \
		int other = swap + 1; \
		if ((other < heap->size) && (heap->cmp_factor * heap->cmp(heap->data[other], heap->data[swap]) > 0)) \
      swap = other; \
		if (heap->cmp_factor * heap->cmp(temp, heap->data[swap]) > 0) break; \
		heap->data[index] = heap->data[swap]; \
		heap->dtors[index] = heap->dtors[swap]; \
	} \
	heap->data[index] = temp; \
	heap->dtors[index] = temp_dtor; \
} \
\
static inline heap_name##_element_t heap_name##_front(heap_name##_t* heap) \
{ \
  return (*heap->data); \
} \
\
static inline bool heap_name##_empty(heap_name##_t* heap) \
{ \
  return (heap->size == 0); \
} \
\
static inline void heap_name##_clear(heap_name##_t* heap) \
{ \
  heap->size = 0; \
  if (heap->capacity > 4) \
  { \
    heap->capacity = 4; \
    heap->data = polymec_realloc(heap->data, sizeof(element) * heap->capacity); \
  } \
} \
\

// Define some heaps.
DEFINE_HEAP(int_heap, int, int_cmp)
DEFINE_HEAP(index_heap, index_t, index_cmp)

#endif