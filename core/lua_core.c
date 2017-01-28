// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/options.h"
#include "core/lua_core.h"
#include "core/partition_mesh.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static int p_new(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 3) || 
      !lua_isnumber(L, 1) || !lua_isnumber(L, 2) || !lua_isnumber(L, 3))
  {
    return luaL_error(L, "Arguments must be x, y, z coordinates.");
  }

  real_t x = (real_t)lua_tonumber(L, 1);
  real_t y = (real_t)lua_tonumber(L, 2);
  real_t z = (real_t)lua_tonumber(L, 3);
  point_t* point = point_new(x, y, z);
  lua_push_point(L, point);
  return 1;
}

static int p_distance(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  point_t* q = lua_to_point(L, 2);
  if (q == NULL)
    return luaL_error(L, "Argument must be a point.");
  lua_pushnumber(L, (double)point_distance(p, q));
  return 1;
}

static lua_module_function point_funcs[] = {
  {"new", p_new},
  {"distance", p_distance},
  {NULL, NULL}
};

static int p_x(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  lua_pushnumber(L, (double)p->x);
  return 1;
}

static int p_set_x(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Point coordinates must be numbers.");
  p->x = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int p_y(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  lua_pushnumber(L, (double)p->y);
  return 1;
}

static int p_set_y(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Point coordinates must be numbers.");
  p->y = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int p_z(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  lua_pushnumber(L, (double)p->z);
  return 1;
}

static int p_set_z(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Point coordinates must be numbers.");
  p->z = (real_t)lua_tonumber(L, 2);
  return 0;
}

static lua_record_field point_fields[] = {
  {"x", p_x, p_set_x},
  {"y", p_y, p_set_y},
  {"z", p_z, p_set_z},
  {NULL, NULL, NULL}
};

static int p_len(lua_State* L)
{
  lua_pushnumber(L, 3.0);
  return 1;
}

static int p_tostring(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  lua_pushfstring(L, "point (%f, %f, %f)", p->x, p->y, p->z);
  return 1;
}

static lua_record_metamethod point_mm[] = {
  {"__len", p_len},
  {"__tostring", p_tostring},
  {NULL, NULL}
};

static int v_new(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 3) || 
      !lua_isnumber(L, 1) || !lua_isnumber(L, 2) || !lua_isnumber(L, 3))
  {
    return luaL_error(L, "Arguments must be x, y, z components.");
  }

  real_t x = (real_t)lua_tonumber(L, 1);
  real_t y = (real_t)lua_tonumber(L, 2);
  real_t z = (real_t)lua_tonumber(L, 3);
  vector_t* vec = vector_new(x, y, z);
  lua_push_vector(L, vec);
  return 1;
}

static int v_dot(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* w = lua_to_vector(L, 2);
  if (w == NULL)
    return luaL_error(L, "Argument must be a vector.");
  lua_pushnumber(L, (double)vector_dot(v, w));
  return 1;
}

static lua_module_function vector_funcs[] = {
  {"new", v_new},
  {"dot", v_dot},
  {NULL, NULL}
};

static int v_x(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  lua_pushnumber(L, (double)v->x);
  return 1;
}

static int v_set_x(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Vector components must be numbers.");
  v->x = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int v_y(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  lua_pushnumber(L, (double)v->y);
  return 1;
}

static int v_set_y(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Vector components must be numbers.");
  v->y = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int v_z(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  lua_pushnumber(L, (double)v->z);
  return 1;
}

static int v_set_z(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Vector components must be numbers.");
  v->z = (real_t)lua_tonumber(L, 2);
  return 0;
}

static lua_record_field vector_fields[] = {
  {"x", v_x, v_set_x},
  {"y", v_y, v_set_y},
  {"z", v_z, v_set_z},
  {NULL, NULL, NULL}
};

static int v_add(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* w = lua_to_vector(L, 2);
  if (v == NULL)
    luaL_error(L, "Argument 1 must be a vector.");
  if (w == NULL)
    luaL_error(L, "Argument 2 must be a vector.");
  vector_t* sum = vector_new(v->x + w->x, v->y + w->y, v->z + w->z);
  lua_push_vector(L, sum);
  return 1;
}

static int v_sub(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* w = lua_to_vector(L, 2);
  if (v == NULL)
    luaL_error(L, "Argument 1 must be a vector.");
  if (w == NULL)
    luaL_error(L, "Argument 2 must be a vector.");
  vector_t* diff = vector_new(v->x - w->x, v->y - w->y, v->z - w->z);
  lua_push_vector(L, diff);
  return 1;
}

static int v_mul(lua_State* L)
{
  if ((!lua_isnumber(L, 1) || !lua_is_vector(L, 2)) &&
      (!lua_is_vector(L, 1) || !lua_isnumber(L, 2)))
    luaL_error(L, "Arguments must be a vector and a number.");
  vector_t* v = lua_to_vector(L, (lua_isnumber(L, 1)) ? 2 : 1);
  real_t c = (real_t)lua_tonumber(L, (lua_isnumber(L, 1)) ? 1 : 2);
  vector_t* v1 = vector_new(c * v->x, c * v->y, c * v->z);
  lua_push_vector(L, v1);
  return 1;
}

static int v_div(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  if (v == NULL)
    luaL_error(L, "Argument 1 must be a vector.");
  if (!lua_isnumber(L, 2))
    luaL_error(L, "Argument 2 must be a number.");
  real_t c = (real_t)lua_tonumber(L, 2);
  vector_t* v1 = vector_new(v->x/c, v->y/c, v->z/c);
  lua_push_vector(L, v1);
  return 1;
}

static int v_unm(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* v1 = vector_new(-1.0 * v->x, -1.0 * v->y, -1.0 * v->z);
  lua_push_vector(L, v1);
  return 1;
}

static int v_len(lua_State* L)
{
  lua_pushnumber(L, 3.0);
  return 1;
}

static int v_tostring(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  lua_pushfstring(L, "vector (%f, %f, %f)", v->x, v->y, v->z);
  return 1;
}

static lua_record_metamethod vector_mm[] = {
  {"__add", v_add},
  {"__sub", v_sub},
  {"__mul", v_mul},
  {"__div", v_div},
  {"__unm", v_unm},
  {"__len", v_len},
  {"__tostring", v_tostring},
  {"dot", v_dot},
  {NULL, NULL}
};

static int t_new(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_istable(L, 1))
  {
    return luaL_error(L, "Argument must be the tensor's shape.");
  }

  int rank = (int)lua_rawlen(L, 1);
  size_t shape[rank];
  for (int i = 1; i <= rank; ++i)
  {
    lua_rawgeti(L, 1, (lua_Integer)i);
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "Shape array must contain only integers.");
    shape[i-1] = (int)lua_tointeger(L, -1);
  }
  tensor_t* t = tensor_new(rank, shape);
  lua_push_tensor(L, t);
  return 1;
}

static lua_module_function tensor_funcs[] = {
  {"new", t_new},
  {NULL, NULL}
};

static int t_rank(lua_State* L)
{
  tensor_t* t = lua_to_tensor(L, 1);
  lua_pushinteger(L, tensor_rank(t));
  return 1;
}

static int t_shape(lua_State* L)
{
  tensor_t* t = lua_to_tensor(L, 1);
  int rank = tensor_rank(t);
  lua_createtable(L, rank, 0);
  size_t* shape = tensor_shape(t);
  for (int i = 1; i <= rank; ++i)
  {
    lua_pushinteger(L, (lua_Integer)i);
    lua_pushinteger(L, (lua_Integer)(shape[i-1]));
    lua_rawset(L, -3);
  }
  return 1;
}

static lua_record_field tensor_fields[] = {
  {"rank", t_rank, NULL},
  {"shape", t_shape, NULL},
  {NULL, NULL, NULL}
};

static void check_tensors_are_same(lua_State* L, tensor_t* t, tensor_t* u)
{
  if (t == NULL)
    luaL_error(L, "Argument 1 must be a tensor.");
  if (u == NULL)
    luaL_error(L, "Argument 2 must be a tensor.");
  int rank = tensor_rank(t);
  if (tensor_rank(u) != rank)
    luaL_error(L, "Arguments 1 and 2 must have matching ranks.");
  size_t* shape = tensor_shape(t);
  size_t* shape1 = tensor_shape(u);
  for (int i = 0; i < rank; ++i)
  {
    if (shape[i] != shape1[i])
      luaL_error(L, "Arguments 1 and 2 must have matching shapes.");
  }
}

static int t_call(lua_State* L)
{
  tensor_t* t = lua_to_tensor(L, 1);
  int rank = tensor_rank(t);
  int num_args = lua_gettop(L);
  if (num_args < (rank + 1))
    return luaL_error(L, "Tensor has %d indices (only %d given).", rank, num_args-1);
  else if (num_args > (rank + 2))
    return luaL_error(L, "Too many arguments. Must use %d indices and possibly one value to assign.", rank);
  size_t offset = 0;
  size_t* shape = tensor_shape(t);
  for (int i = 1; i <= rank; ++i)
  {
    if (!lua_isinteger(L, 1+i))
      return luaL_error(L, "Tensor index %d is not an integer.", i);
    int index = (int)(lua_tointeger(L, 1+i));
    if ((index < 1) || (index > shape[i-1]))
      return luaL_error(L, "Tensor index %d is out of range: %d.", i, index);
    if (i == rank)
      offset += (index-1);
    else
      offset += (index-1) * shape[i-1];
  }
  real_t* t_data = tensor_data(t);
  if (num_args == (rank + 2)) // set
  {
    if (!lua_isnumber(L, num_args))
      return luaL_error(L, "Assigned tensor value must be a number.");
    t_data[offset] = (real_t)lua_tonumber(L, num_args);
    return 0;
  }
  else // get
  {
    lua_pushnumber(L, t_data[offset]);
    return 1;
  }
}

static int t_add(lua_State* L)
{
  tensor_t* t = lua_to_tensor(L, 1);
  tensor_t* u = lua_to_tensor(L, 2);
  check_tensors_are_same(L, t, u);
  tensor_t* sum = tensor_clone(t);
  real_t* u_data = tensor_data(u);
  real_t* sum_data = tensor_data(sum);
  size_t size = tensor_size(t);
  for (size_t i = 0; i < size; ++i)
    sum_data[i] += u_data[i];
  lua_push_tensor(L, sum);
  return 1;
}

static int t_sub(lua_State* L)
{
  tensor_t* t = lua_to_tensor(L, 1);
  tensor_t* u = lua_to_tensor(L, 2);
  check_tensors_are_same(L, t, u);
  tensor_t* diff = tensor_clone(t);
  real_t* u_data = tensor_data(u);
  real_t* diff_data = tensor_data(diff);
  size_t size = tensor_size(t);
  for (size_t i = 0; i < size; ++i)
    diff_data[i] -= u_data[i];
  lua_push_tensor(L, diff);
  return 1;
}

static int t_mul(lua_State* L)
{
  if ((!lua_isnumber(L, 1) || !lua_is_tensor(L, 2)) &&
      (!lua_is_tensor(L, 1) || !lua_isnumber(L, 2)))
    luaL_error(L, "Arguments must be a tensor and a number.");
  tensor_t* t = lua_to_tensor(L, (lua_isnumber(L, 1)) ? 2 : 1);
  real_t c = (real_t)lua_tonumber(L, (lua_isnumber(L, 1)) ? 1 : 2);
  tensor_t* t1 = tensor_clone(t);
  real_t* t1_data = tensor_data(t1);
  size_t size = tensor_size(t);
  for (size_t i = 0; i < size; ++i)
    t1_data[i] *= c;
  lua_push_tensor(L, t1);
  return 1;
}

static int t_div(lua_State* L)
{
  tensor_t* t = lua_to_tensor(L, 1);
  if (t == NULL)
    luaL_error(L, "Argument 1 must be a tensor.");
  if (!lua_isnumber(L, 2))
    luaL_error(L, "Argument 2 must be a number.");
  real_t c = (real_t)lua_tonumber(L, 2);
  tensor_t* t1 = tensor_clone(t);
  real_t* t1_data = tensor_data(t1);
  size_t size = tensor_size(t);
  for (size_t i = 0; i < size; ++i)
    t1_data[i] /= c;
  lua_push_tensor(L, t1);
  return 1;
}

static int t_unm(lua_State* L)
{
  tensor_t* t = lua_to_tensor(L, 1);
  tensor_t* t1 = tensor_clone(t);
  real_t* t1_data = tensor_data(t1);
  size_t size = tensor_size(t);
  for (size_t i = 0; i < size; ++i)
    t1_data[i] = -t1_data[i];
  lua_push_tensor(L, t1);
  return 1;
}

static int t_len(lua_State* L)
{
  tensor_t* t = lua_to_tensor(L, 1);
  lua_pushinteger(L, tensor_rank(t));
  return 1;
}

static int t_tostring(lua_State* L)
{
  tensor_t* t = lua_to_tensor(L, 1);
  int rank = tensor_rank(t);
  if (rank == 0)
    lua_pushfstring(L, "tensor (rank %d, value = %f)", tensor_rank(t), tensor_data(t)[0]);
  else
    lua_pushfstring(L, "tensor (rank %d)", tensor_rank(t));
  return 1;
}

static lua_record_metamethod tensor_mm[] = {
  {"__call", t_call},
  {"__add", t_add},
  {"__sub", t_sub},
  {"__mul", t_mul},
  {"__div", t_div},
  {"__unm", t_unm},
  {"__len", t_len},
  {"__tostring", t_tostring},
  {NULL, NULL}
};

static int bb_new(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if (num_args != 1)
  {
    if (!lua_istable(L, 1))
      return luaL_error(L, "Argument must be a table containing x1, x2, y1, y2, z1, z2 values.");
  }

  // Look for x1, x2, y1, y2, z1, z2 in the table.
  bbox_t* bbox = bbox_new(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  const char* entries[] = {"x1", "x2", "y1", "y2", "z1", "z2"};
  for (int i = 0; i < 6; ++i)
  {
    lua_pushstring(L, entries[i]);
    lua_gettable(L, 1); // Reads name from top, replaces with bounds[name].
    if (!lua_isnumber(L, -1))
    {
      return luaL_error(L, "Invalid entry for '%s'.\n"
                        "x1, x2, y1, y2, z1, z2, must all be numbers.", entries[i]);
    }
    switch(i)
    {
      case 0: bbox->x1 = (real_t)lua_tonumber(L, -1);
              break;
      case 1: bbox->x2 = (real_t)lua_tonumber(L, -1);
              break;
      case 2: bbox->y1 = (real_t)lua_tonumber(L, -1);
              break;
      case 3: bbox->y2 = (real_t)lua_tonumber(L, -1);
              break;
      case 4: bbox->z1 = (real_t)lua_tonumber(L, -1);
              break;
      case 5: bbox->z2 = (real_t)lua_tonumber(L, -1);
              break;
      default: break;
    }
    lua_pop(L, 1); 
  }

  // Push the bounding box onto the stack.
  lua_push_bbox(L, bbox);
  return 1;
}

static int bb_contains(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  point_t* p = lua_to_point(L, 2);
  if (p == NULL)
    return luaL_error(L, "Argument must be a point.");
  lua_pushboolean(L, bbox_contains(b, p));
  return 1;
}

static lua_module_function bbox_funcs[] = {
  {"new", bb_new},
  {NULL, NULL}
};

static int bb_tostring(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  lua_pushfstring(L, "bbox (x1 = %f, x2 = %f, y1 = %f, y2 = %f, z1 = %f, z2 = %f)", 
                  b->x1, b->x2, b->y1, b->y2, b->z1, b->z2);
  return 1;
}

static lua_class_method bbox_methods[] = {
  {"contains", bb_contains},
  {"__tostring", bb_tostring},
  {NULL, NULL}
};

static int sp_constant(lua_State* L)
{
  // Check the argument.
  int num_args = lua_gettop(L);
  real_t val[num_args];
  for (int i = 0; i < num_args; ++i)
  {
    if (!lua_isnumber(L, i))
      return luaL_error(L, "Argument %d must be a number.", i);
    val[i] = (real_t)lua_tonumber(L, i);
  }
  sp_func_t* f = constant_sp_func_new(val, num_args);
  lua_push_sp_func(L, f);
  return 1;
}

static lua_module_function sp_funcs[] = {
  {"constant", sp_constant},
  {NULL, NULL}
};

static int sp_len(lua_State* L)
{
  sp_func_t* f = lua_to_sp_func(L, 1);
  lua_pushinteger(L, sp_func_num_comp(f));
  return 1;
}

static int sp_call(lua_State* L)
{
  sp_func_t* f = lua_to_sp_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument must be a point.");
  point_t* x = lua_to_point(L, 2);
  int nc = sp_func_num_comp(f);
  real_t val[nc];
  sp_func_eval(f, x, val);
  for (int i = 0; i < nc; ++i)
    lua_pushnumber(L, val[i]);
  return nc;
}

static lua_class_method sp_methods[] = {
  {"__len", sp_len},
  {"__call", sp_call},
  {NULL, NULL}
};

static int st_constant(lua_State* L)
{
  // Check the argument.
  int num_args = lua_gettop(L);
  real_t val[num_args];
  for (int i = 0; i < num_args; ++i)
  {
    if (!lua_isnumber(L, i))
      return luaL_error(L, "Argument %d must be a number.", i);
    val[i] = (real_t)lua_tonumber(L, i);
  }
  st_func_t* f = constant_st_func_new(val, num_args);
  lua_push_st_func(L, f);
  return 1;
}

static lua_module_function st_funcs[] = {
  {"constant", st_constant},
  {NULL, NULL}
};

static int st_len(lua_State* L)
{
  st_func_t* f = lua_to_st_func(L, 1);
  lua_pushinteger(L, st_func_num_comp(f));
  return 1;
}

static int st_call(lua_State* L)
{
  st_func_t* f = lua_to_st_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "First argument must be a point.");
  if (!lua_is_point(L, 3))
    return luaL_error(L, "Second argument must be a time.");
  point_t* x = lua_to_point(L, 2);
  real_t t = (real_t)lua_tonumber(L, 3);
  int nc = st_func_num_comp(f);
  real_t val[nc];
  st_func_eval(f, x, t, val);
  for (int i = 0; i < nc; ++i)
    lua_pushnumber(L, val[i]);
  return nc;
}

static lua_class_method st_methods[] = {
  {"__len", st_len},
  {"__call", st_call},
  {NULL, NULL}
};

static int mesh_repartition(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_is_mesh(L, 1))
    return luaL_error(L, "Argument must be mesh.");
  mesh_t* mesh = lua_to_mesh(L, 1);
  real_t imbalance_tol = 0.05;

  // Bug out if there's only one process.
  int nprocs;
  MPI_Comm_size(mesh->comm, &nprocs);
  if (nprocs == 1)
    return 0;

  // Make sure there are enough cells for our processes.
  index_t local_num_cells = mesh->num_cells, global_num_cells;
  MPI_Allreduce(&local_num_cells, &global_num_cells, 1, MPI_INDEX_T, MPI_SUM, mesh->comm);
  if (global_num_cells < nprocs)
    return luaL_error(L, "Insufficient number of cells (%zd) for number of processes (%d).", global_num_cells, nprocs);

  // Perform the repartitioning and toss the migrator.
  migrator_t* m = repartition_mesh(&mesh, NULL, imbalance_tol);
  m = NULL;

  return 0;
}

static lua_module_function mesh_funcs[] = {
  {"repartition", mesh_repartition},
  {NULL, NULL}
};

static int mesh_num_cells(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushinteger(L, m->num_cells);
  return 1;
}

static int mesh_num_ghost_cells(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushinteger(L, m->num_ghost_cells);
  return 1;
}

static int mesh_num_faces(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushinteger(L, m->num_faces);
  return 1;
}

static int mesh_num_edges(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushinteger(L, m->num_edges);
  return 1;
}

static int mesh_num_nodes(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushinteger(L, m->num_nodes);
  return 1;
}

static lua_record_field mesh_fields[] = {
  {"num_cells", mesh_num_cells, NULL},
  {"num_ghost_cells", mesh_num_ghost_cells, NULL},
  {"num_faces", mesh_num_faces, NULL},
  {"num_edges", mesh_num_edges, NULL},
  {"num_nodes", mesh_num_nodes, NULL},
  {NULL, NULL, NULL}
};

static int mesh_tostring(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushfstring(L, "mesh (%d cells, %d faces, %d nodes)", 
                  m->num_cells, m->num_faces, m->num_nodes);
  return 1;
}

static lua_record_metamethod mesh_mm[] = {
  {"__tostring", mesh_tostring},
  {NULL, NULL}
};

static int pc_repartition(lua_State* L)
{
  luaL_error(L, "can't repartition point clouds just yet!");
  return 0;
}

static lua_module_function pc_funcs[] = {
  {"repartition", pc_repartition},
  {NULL, NULL}
};

static int pc_num_points(lua_State* L)
{
  point_cloud_t* pc = lua_to_point_cloud(L, 1);
  lua_pushinteger(L, pc->num_points);
  return 1;
}

static int pc_num_ghosts(lua_State* L)
{
  point_cloud_t* pc = lua_to_point_cloud(L, 1);
  lua_pushinteger(L, pc->num_ghosts);
  return 1;
}

static lua_record_field pc_fields[] = {
  {"num_points", pc_num_points, NULL},
  {"num_ghosts", pc_num_ghosts, NULL},
  {NULL, NULL, NULL}
};

static int pc_tostring(lua_State* L)
{
  point_cloud_t* pc = lua_to_point_cloud(L, 1);
  lua_pushfstring(L, "point cloud (%d points)", pc->num_points);
  return 1;
}

static lua_record_metamethod pc_mm[] = {
  {"__len", pc_num_points},
  {"__tostring", pc_tostring},
  {NULL, NULL}
};

static int silo_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with prefix, dir entries.");

  lua_getfield(L, -2, "prefix");
  if (lua_isnil(L, -1))
    return luaL_error(L, "prefix must be specified.");
  const char* prefix = lua_tostring(L, -1);
  lua_pop(L, 1);

  lua_getfield(L, -2, "dir");
  if (lua_isnil(L, -1))
    return luaL_error(L, "dir must be specified.");
  const char* dir = lua_tostring(L, -1);
  lua_pop(L, 1);

  int num_files = 1;
  lua_getfield(L, -2, "num_files");
  if (!lua_isnil(L, -1))
    num_files = (int)(lua_tonumber(L, -1));
  lua_pop(L, 1);

  int step = 0;
  lua_getfield(L, -2, "step");
  if (!lua_isnil(L, -1))
    step = (int)(lua_tonumber(L, -1));
  lua_pop(L, 1);

  real_t time = 0.0;
  lua_getfield(L, -2, "time");
  if (!lua_isnil(L, -1))
    time = (real_t)(lua_tonumber(L, -1));
  lua_pop(L, 1);

  silo_file_t* s = silo_file_new(MPI_COMM_WORLD, prefix, dir, num_files, 0, step, time);
  lua_push_silo_file(L, s);
  return 1;
}

static int silo_open(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with prefix, dir, step entries.");

  lua_getfield(L, -2, "prefix");
  if (lua_isnil(L, -1))
    return luaL_error(L, "prefix must be specified.");
  const char* prefix = lua_tostring(L, -1);
  lua_pop(L, 1);

  lua_getfield(L, -2, "dir");
  if (lua_isnil(L, -1))
    return luaL_error(L, "dir must be specified.");
  const char* dir = lua_tostring(L, -1);
  lua_pop(L, 1);

  int num_files = 1;
  lua_getfield(L, -2, "num_files");
  if (!lua_isnil(L, -1))
    num_files = (int)(lua_tonumber(L, -1));
  lua_pop(L, 1);

  int step = 0;
  lua_getfield(L, -2, "step");
  if (lua_isnil(L, -1))
    return luaL_error(L, "step must be specified.");
  step = (int)(lua_tonumber(L, -1));
  lua_pop(L, 1);

  real_t time;
  silo_file_t* s = silo_file_open(MPI_COMM_WORLD, prefix, dir, 0, step, &time);
  lua_push_silo_file(L, s);
  return 1;
}

static int silo_close(lua_State* L)
{
  silo_file_t* s = lua_to_silo_file(L, 1);
  silo_file_close(s);
  return 0;
}

static lua_module_function silo_funcs[] = {
  {"new", silo_new},
  {"open", silo_open},
  {NULL, NULL}
};

static lua_class_method silo_methods[] = {
  {"close", silo_close},
  {NULL, NULL}
};

static void register_options(lua_State* L)
{
  // Create a new table and fill it with our named command line values.
  lua_newtable(L);
  options_t* opts = options_argv();
  int pos = 0;
  const char* opt_name;
  const char* opt_val;
  while (options_next_value(opts, &pos, &opt_name, &opt_val))
  {
    lua_pushstring(L, opt_val);
    lua_setfield(L, -2, opt_name);
  }
  lua_setglobal(L, "options");
}

//------------------------------------------------------------------------
//                                API 
//------------------------------------------------------------------------

int lua_register_core_modules(lua_State* L)
{
  // Core types.
  lua_register_record_type(L, "point", point_funcs, point_fields, point_mm);
  lua_register_record_type(L, "vector", vector_funcs, vector_fields, vector_mm);
  lua_register_record_type(L, "tensor", tensor_funcs, tensor_fields, tensor_mm);
  lua_register_class(L, "bbox", bbox_funcs, bbox_methods);
  lua_register_class(L, "sp_func", sp_funcs, sp_methods);
  lua_register_class(L, "st_func", st_funcs, st_methods);
  lua_register_record_type(L, "mesh", mesh_funcs, mesh_fields, mesh_mm);
  lua_register_record_type(L, "point_cloud", pc_funcs, pc_fields, pc_mm);
  lua_register_class(L, "silo_file", silo_funcs, silo_methods);

  // Register the options table.
  register_options(L);

  return 0;
}

void lua_push_point(lua_State* L, point_t* p)
{
  lua_push_record(L, "point", p, NULL);
}

bool lua_is_point(lua_State* L, int index)
{
  return lua_is_record(L, index, "point");
}

point_t* lua_to_point(lua_State* L, int index)
{
  return (point_t*)lua_to_record(L, index, "point");
}

bool lua_is_point_list(lua_State* L, int index)
{
  if (lua_istable(L, index)) 
  {
    size_t len = lua_rawlen(L, index);
    if (len == 0) 
      return false;
    bool is_point_list = true;
    for (size_t i = 1; i <= len; ++i)
    {
      lua_rawgeti(L, index, (lua_Integer)i);
      bool is_point = lua_is_point(L, -1);
      if (!is_point)
      {
        bool is_3_tuple = (lua_istable(L, -1) && 
                           (lua_rawlen(L, -1) == 3));
        if (!is_3_tuple)
          is_point_list = false;
        else
        {
          for (size_t j = 1; j <= 3; ++j)
          {
            lua_rawgeti(L, index, (lua_Integer)j);
            if (!lua_isnumber(L, -1))
              is_point_list = false;
            lua_pop(L, 1);
            if (!is_point_list)
              break;
          }
        }
        if (!is_point_list)
          break;
      }
    }
    return is_point_list;
  }
  else 
    return false;
}

bool lua_is_canonical_point_list(lua_State* L, int index)
{
  if (lua_istable(L, index)) 
  {
    size_t len = lua_rawlen(L, index);
    if (len == 0) 
      return false;
    bool is_point_list = true;
    for (size_t i = 1; i <= len; ++i)
    {
      lua_rawgeti(L, index, (lua_Integer)i);
      bool is_point = lua_is_point(L, -1);
      if (!is_point)
      {
        is_point_list = false;
        break;
      }
    }
    return is_point_list;
  }
  else 
    return false;
}

void lua_canonicalize_point_list(lua_State* L, int index)
{
  if (lua_is_point_list(L, index))
  {
    size_t len = lua_rawlen(L, index);
    for (size_t i = 1; i <= len; ++i)
    {
      lua_rawgeti(L, index, (lua_Integer)i);
      if (!lua_is_point(L, -1))
      {
        real_t pj[3];
        for (int j = 1; j <= 3; ++j)
        {
          lua_rawgeti(L, -1, (lua_Integer)j);
          pj[j-1] = (real_t)lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
        point_t* p = point_new(pj[0], pj[1], pj[2]);
        lua_push_point(L, p);
        lua_rawseti(L, -2, i); 
      }
    }
  }
}

void lua_export_point_list(lua_State* L, int index, real_t* array)
{
  if (lua_is_canonical_point_list(L, index))
  {
    size_t len = lua_rawlen(L, index);
    for (size_t i = 1; i <= len; ++i)
    {
      lua_rawgeti(L, index, (lua_Integer)i);
      point_t* pi = lua_to_point(L, -1);
      array[3*(i-1)  ] = pi->x;
      array[3*(i-1)+1] = pi->y;
      array[3*(i-1)+2] = pi->z;
    }
  }
}

void lua_import_point_list(lua_State* L, int index, real_t* array)
{
  if (lua_is_canonical_point_list(L, index))
  {
    size_t len = lua_rawlen(L, index);
    for (size_t i = 1; i <= len; ++i)
    {
      lua_rawgeti(L, index, (lua_Integer)i);
      point_t* pi = lua_to_point(L, -1);
      pi->x = array[3*(i-1)];
      pi->y = array[3*(i-1)+1];
      pi->z = array[3*(i-1)+2];
    }
  }
}

void lua_push_vector(lua_State* L, vector_t* v)
{
  lua_push_record(L, "vector", v, NULL);
}

bool lua_is_vector(lua_State* L, int index)
{
  return lua_is_record(L, index, "vector");
}

vector_t* lua_to_vector(lua_State* L, int index)
{
  return (vector_t*)lua_to_record(L, index, "vector");
}

bool lua_is_vector_list(lua_State* L, int index)
{
  if (lua_istable(L, index)) 
  {
    size_t len = lua_rawlen(L, index);
    if (len == 0) 
      return false;
    bool is_vector_list = true;
    for (size_t i = 1; i <= len; ++i)
    {
      lua_rawgeti(L, index, (lua_Integer)i);
      bool is_vector = lua_is_vector(L, -1);
      if (!is_vector)
      {
        bool is_3_tuple = (lua_istable(L, -1) && 
                           (lua_rawlen(L, -1) == 3));
        if (!is_3_tuple)
          is_vector_list = false;
        else
        {
          for (size_t j = 1; j <= 3; ++j)
          {
            lua_rawgeti(L, index, (lua_Integer)j);
            if (!lua_isnumber(L, -1))
              is_vector_list = false;
            lua_pop(L, 1);
            if (!is_vector_list)
              break;
          }
        }
        if (!is_vector_list)
          break;
      }
    }
    return is_vector_list;
  }
  else 
    return false;
}

bool lua_is_canonical_vector_list(lua_State* L, int index)
{
  if (lua_istable(L, index)) 
  {
    size_t len = lua_rawlen(L, index);
    if (len == 0) 
      return false;
    bool is_vector_list = true;
    for (size_t i = 1; i <= len; ++i)
    {
      lua_rawgeti(L, index, (lua_Integer)i);
      bool is_vector = lua_is_vector(L, -1);
      if (!is_vector)
      {
        is_vector_list = false;
        break;
      }
    }
    return is_vector_list;
  }
  else 
    return false;
}

void lua_canonicalize_vector_list(lua_State* L, int index)
{
  if (lua_is_vector_list(L, index))
  {
    size_t len = lua_rawlen(L, index);
    for (size_t i = 1; i <= len; ++i)
    {
      lua_rawgeti(L, index, (lua_Integer)i);
      if (!lua_is_vector(L, -1))
      {
        real_t vj[3];
        for (int j = 1; j <= 3; ++j)
        {
          lua_rawgeti(L, -1, (lua_Integer)j);
          vj[j-1] = (real_t)lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
        vector_t* v = vector_new(vj[0], vj[1], vj[2]);
        lua_push_vector(L, v);
        lua_rawseti(L, -2, i); 
      }
    }
  }
}

void lua_export_vector_list(lua_State* L, int index, real_t* array)
{
  if (lua_is_canonical_vector_list(L, index))
  {
    size_t len = lua_rawlen(L, index);
    for (size_t i = 1; i <= len; ++i)
    {
      lua_rawgeti(L, index, (lua_Integer)i);
      vector_t* vi = lua_to_vector(L, -1);
      array[3*(i-1)  ] = vi->x;
      array[3*(i-1)+1] = vi->y;
      array[3*(i-1)+2] = vi->z;
    }
  }
}

void lua_import_vector_list(lua_State* L, int index, real_t* array)
{
  if (lua_is_canonical_vector_list(L, index))
  {
    size_t len = lua_rawlen(L, index);
    for (size_t i = 1; i <= len; ++i)
    {
      lua_rawgeti(L, index, (lua_Integer)i);
      vector_t* vi = lua_to_vector(L, -1);
      vi->x = array[3*(i-1)];
      vi->y = array[3*(i-1)+1];
      vi->z = array[3*(i-1)+2];
    }
  }
}

void lua_push_bbox(lua_State* L, bbox_t* b)
{
  lua_push_object(L, "bbox", b, NULL);
}

bool lua_is_bbox(lua_State* L, int index)
{
  return lua_is_object(L, index, "bbox");
}

bbox_t* lua_to_bbox(lua_State* L, int index)
{
  return (bbox_t*)lua_to_object(L, index, "bbox");
}

void lua_push_sp_func(lua_State* L, sp_func_t* f)
{
  lua_push_object(L, "sp_func", f, NULL);
}

bool lua_is_sp_func(lua_State* L, int index)
{
  return lua_is_object(L, index, "sp_func");
}

sp_func_t* lua_to_sp_func(lua_State* L, int index)
{
  return (sp_func_t*)lua_to_object(L, index, "sp_func");
}

void lua_push_st_func(lua_State* L, st_func_t* f)
{
  lua_push_object(L, "st_func", f, NULL);
}

bool lua_is_st_func(lua_State* L, int index)
{
  return lua_is_object(L, index, "st_func");
}

st_func_t* lua_to_st_func(lua_State* L, int index)
{
  return (st_func_t*)lua_to_object(L, index, "st_func");
}

void lua_push_mesh(lua_State* L, mesh_t* m)
{
  lua_push_object(L, "mesh", m, DTOR(mesh_free));
}

bool lua_is_mesh(lua_State* L, int index)
{
  return lua_is_object(L, index, "mesh");
}

mesh_t* lua_to_mesh(lua_State* L, int index)
{
  return (mesh_t*)lua_to_object(L, index, "mesh");
}

void lua_push_point_cloud(lua_State* L, point_cloud_t* c)
{
  lua_push_object(L, "point_cloud", c, DTOR(point_cloud_free));
}

bool lua_is_point_cloud(lua_State* L, int index)
{
  return lua_is_object(L, index, "point_cloud");
}

point_cloud_t* lua_to_point_cloud(lua_State* L, int index)
{
  return (point_cloud_t*)lua_to_object(L, index, "point_cloud");
}

void lua_push_silo_file(lua_State* L, silo_file_t* s)
{
  lua_push_object(L, "silo_file", s, DTOR(silo_file_close));
}

bool lua_is_silo_file(lua_State* L, int index)
{
  return lua_is_object(L, index, "silo_file");
}

silo_file_t* lua_to_silo_file(lua_State* L, int index)
{
  return (silo_file_t*)lua_to_object(L, index, "silo_file");
}

void lua_push_tensor(lua_State* L, tensor_t* t)
{
  lua_push_record(L, "tensor", t, DTOR(tensor_free));
}

bool lua_is_tensor(lua_State* L, int index)
{
  return lua_is_record(L, index, "tensor");
}

tensor_t* lua_to_tensor(lua_State* L, int index)
{
  return (tensor_t*)lua_to_record(L, index, "tensor");
}

