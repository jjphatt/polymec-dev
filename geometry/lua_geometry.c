// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/lua_core.h"
#include "geometry/lua_geometry.h"
#include "geometry/plane_sd_func.h"
#include "geometry/sphere_sd_func.h"
#include "geometry/cylinder_sd_func.h"
#include "geometry/union_sd_func.h"
#include "geometry/intersection_sd_func.h"
#include "geometry/difference_sd_func.h"

#include "geometry/create_uniform_mesh.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// This serves as a context pointer that houses Lua objects.
typedef struct
{
  lua_State* L;
} lua_obj_t;

static real_t sd_value(void* context, point_t* x)
{
  // Fetch our Lua table from the registry.
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);

  // Call the value method.
  lua_getfield(lo->L, -1, "value");
  lua_pushvalue(lo->L, -2);
  lua_push_point(lo->L, x);
  lua_call(lo->L, 2, 1);

  if (!lua_isnumber(lo->L, -1))
    return luaL_error(lo->L, "value method did not return a number.");
  return lua_to_real(lo->L, -1);
}

static void sd_eval_grad(void* context, point_t* x, vector_t* grad)
{
  // Fetch our Lua table from the registry.
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);

  // Call the grad method.
  lua_getfield(lo->L, -1, "grad");
  lua_pushvalue(lo->L, -2);
  lua_push_point(lo->L, x);
  lua_call(lo->L, 2, 1);

  if (!lua_is_vector(lo->L, -1))
    luaL_error(lo->L, "grad method did not return a vector.");
  *grad = *lua_to_vector(lo->L, -1);
}

static void sd_dtor(void* context)
{
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);
  lua_pushnil(lo->L);
  lua_settable(lo->L, LUA_REGISTRYINDEX);
  polymec_free(context);
}

static int sd_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with a name and methods.");

  // Create a context that houses this object.
  lua_newtable(L);
  int obj_index = num_args + 1;

  sd_func_vtable vtable = {.dtor = sd_dtor};

  lua_getfield(L, 1, "name");
  if (!lua_isstring(L, -1))
    return luaL_error(L, "Name must be a string.");
  const char* name = lua_tostring(L, -1);

  lua_getfield(L, 1, "value");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "value must be a method.");
  lua_setfield(L, obj_index, "value");
  vtable.value = sd_value;

  lua_getfield(L, 1, "grad");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "grad must be a method.");
  lua_setfield(L, obj_index, "grad");
  vtable.eval_grad = sd_eval_grad;

  // Allocate a context pointer and stuff our object into the registry.
  lua_obj_t* lo = polymec_malloc(sizeof(lua_obj_t));
  lo->L = L;

  // Store the table representing our object in the registry, with 
  // context as a key.
  lua_pushvalue(L, obj_index);
  lua_rawsetp(L, LUA_REGISTRYINDEX, lo);

  // Set up the function.
  sd_func_t* f = sd_func_new(name, lo, vtable);
  lua_push_sd_func(L, f);

  return 1;
}

static int sd_from_sp_funcs(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 3)
    return luaL_error(L, "Arguments must be a name, a signed distance function and its gradient.");
  if (!lua_isstring(L, 1))
    return luaL_error(L, "Argument 1 must be a name.");
  if (!lua_is_sp_func(L, 2))
    return luaL_error(L, "Argument 2 must be a signed distance function.");
  if (!lua_is_sp_func(L, 3))
    return luaL_error(L, "Argument 3 must be the gradient of argument 2.");
  const char* name = lua_tostring(L, 1);
  sp_func_t* D = lua_to_sp_func(L, 2);
  sp_func_t* G = lua_to_sp_func(L, 3);
  lua_push_sd_func(L, sd_func_from_sp_funcs(name, D, G));
  return 1;
}

static int sd_plane(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) && !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with normal and point entries.");

  lua_getfield(L, 1, "normal");
  if (!lua_is_vector(L, -1))
    return luaL_error(L, "normal must be a vector.");
  vector_t* n = lua_to_vector(L, -1);

  lua_getfield(L, 1, "point");
  if (!lua_is_point(L, -1))
    return luaL_error(L, "point must be a point.");
  point_t* x = lua_to_point(L, -1);

  sd_func_t* f = plane_sd_func_new(n, x);
  lua_push_sd_func(L, f);
  return 1;
}

static int sd_sphere(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) && !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with point, radius, and orientation entries.");

  lua_getfield(L, 1, "point");
  if (!lua_is_point(L, -1))
    return luaL_error(L, "point must be a point.");
  point_t* x = lua_to_point(L, -1);

  lua_getfield(L, 1, "radius");
  if (!lua_isnumber(L, -1))
    return luaL_error(L, "radius must be a positive number.");
  real_t r = lua_to_real(L, -1);
  if (r <= 0.0)
    return luaL_error(L, "radius must be positive.");

  normal_orient_t orient;
  lua_getfield(L, 1, "orientation");
  if (!lua_isstring(L, -1))
    return luaL_error(L, "orientation must be 'inward' or 'outward'.");
  const char* orientation = lua_tostring(L, -1);
  if (strcasecmp(orientation, "inward") == 0)
    orient = INWARD_NORMAL;
  else if (strcasecmp(orientation, "outward") == 0)
    orient = OUTWARD_NORMAL;
  else 
    return luaL_error(L, "invalid orientation: %s (must be 'inward' or 'outward').", orientation);

  sd_func_t* f = sphere_sd_func_new(x, r, orient);
  lua_push_sd_func(L, f);
  return 1;
}

static int sd_cylinder(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) && !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with point, radius, and orientation entries.");

  lua_getfield(L, 1, "point");
  if (!lua_is_point(L, -1))
    return luaL_error(L, "point must be a point.");
  point_t* x = lua_to_point(L, -1);

  lua_getfield(L, 1, "radius");
  if (!lua_isnumber(L, -1))
    return luaL_error(L, "radius must be a positive number.");
  real_t r = lua_to_real(L, -1);
  if (r <= 0.0)
    return luaL_error(L, "radius must be positive.");

  normal_orient_t orient;
  lua_getfield(L, 1, "orientation");
  if (!lua_isstring(L, -1))
    return luaL_error(L, "orientation must be 'inward' or 'outward'.");
  const char* orientation = lua_tostring(L, -1);
  if (strcasecmp(orientation, "inward") == 0)
    orient = INWARD_NORMAL;
  else if (strcasecmp(orientation, "outward") == 0)
    orient = OUTWARD_NORMAL;
  else 
    return luaL_error(L, "invalid orientation: %s (must be 'inward' or 'outward').", orientation);

  sd_func_t* f = cylinder_sd_func_new(x, r, orient);
  lua_push_sd_func(L, f);
  return 1;
}

static int sd_union(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) && !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table of sd_funcs.");

  int num_surfaces = (int)lua_rawlen(L, 1);
  sd_func_t* surfaces[num_surfaces];

  for (int i = 0; i < num_surfaces; ++i)
  {
    lua_rawgeti(L, 1, i+1);
    if (!lua_is_sd_func(L, -1))
      return luaL_error(L, "Item %d in table is not an sd_func.", i+1);
    surfaces[i] = lua_to_sd_func(L, -1);
  }

  sd_func_t* f = union_sd_func_new(surfaces, num_surfaces);
  lua_push_sd_func(L, f);
  return 1;
}

static int sd_intersection(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) && !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table of sd_funcs.");

  int num_surfaces = (int)lua_rawlen(L, 1);
  sd_func_t* surfaces[num_surfaces];

  for (int i = 0; i < num_surfaces; ++i)
  {
    lua_rawgeti(L, 1, i+1);
    if (!lua_is_sd_func(L, -1))
      return luaL_error(L, "Item %d in table is not an sd_func.", i+1);
    surfaces[i] = lua_to_sd_func(L, -1);
  }

  sd_func_t* f = intersection_sd_func_new(surfaces, num_surfaces);
  lua_push_sd_func(L, f);
  return 1;
}

static int sd_difference(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2)
    return luaL_error(L, "Arguments must be a pair of sd_funcs.");
  if (!lua_is_sd_func(L, 1))
    return luaL_error(L, "Argument 1 must be an sd_func.");
  if (!lua_is_sd_func(L, 2))
    return luaL_error(L, "Argument 2 must be an sd_func.");

  sd_func_t* f = difference_sd_func_new(lua_to_sd_func(L, 1), 
                                        lua_to_sd_func(L, 2));
  lua_push_sd_func(L, f);
  return 1;
}

static lua_module_function sd_funcs[] = {
  {"new", sd_new, "sd_func.new{name = NAME, value = F, grad = gradF} -> new signed distance function."},
  {"from_sp_funcs", sd_from_sp_funcs, "sd_func.from_sp_funcs(name, D, G) -> new signed distance function using D for distance and G for gradients."},
  {"plane", sd_plane, "sd_func.plane{normal = n, point = x} -> new plane signed distance function."},
  {"sphere", sd_sphere, "sd_func.sphere{point = x, radius = r, orientation = 'inward'|'outward'} -> new sphere signed distance function."},
  {"cylinder", sd_cylinder, "sd_func.cylinder{point = x, radius = r, orientation = 'inward'|'outward'} -> new (infinite) cylinder signed distance function."},
  {"union", sd_union, "sd_func.union(funcs) -> new signed distance function (union of funcs)."},
  {"intersection", sd_intersection, "sd_func.intersection(funcs) -> new signed distance function (intersection of funcs)."},
  {"difference", sd_difference, "sd_func.difference(f, g) -> new signed distance function (f - g)."},
  {NULL, NULL, NULL}
};

static int sd_rename(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a string.");
  sd_func_rename(f, lua_tostring(L, 2));
  return 0;
}

static int sd_grad(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument 1 must be a point.");
  point_t* x = lua_to_point(L, 2);
  vector_t* grad = vector_new(0.0, 0.0, 0.0);
  sd_func_eval_grad(f, x, grad);
  lua_push_vector(L, grad);
  return 1;
}

static int sd_call(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument must be a point.");
  point_t* x = lua_to_point(L, 2);
  lua_pushnumber(L, sd_func_value(f, x));
  return 1;
}

static int sd_project(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument must be a point.");
  point_t* x = lua_to_point(L, 2);
  point_t* proj_x = point_new(0.0, 0.0, 0.0);
  sd_func_project(f, x, proj_x);
  lua_push_point(L, proj_x);
  return 1;
}

static int sd_tostring(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  lua_pushfstring(L, "sd_func '%s'", sd_func_name(f));
  return 1;
}

static lua_class_method sd_methods[] = {
  {"rename", sd_rename, "f:rename(name) -> Renames f to the given name."},
  {"grad", sd_grad, "f:grad(x) -> Returns the gradient of f at x."},
  {"project", sd_project, "f:project(x) -> Returns the projection of f to x."},
  {"__call", sd_call, NULL},
  {"__tostring", sd_tostring, NULL},
  {NULL, NULL, NULL}
};

static real_t sdt_value(void* context, point_t* x, real_t t)
{
  // Fetch our Lua table from the registry.
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);

  // Call the value method.
  lua_getfield(lo->L, -1, "value");
  lua_pushvalue(lo->L, -2);
  lua_push_point(lo->L, x);
  lua_push_real(lo->L, t);
  lua_call(lo->L, 3, 1);

  if (!lua_isnumber(lo->L, -1))
    return luaL_error(lo->L, "value method did not return a number.");
  return lua_to_real(lo->L, -1);
}

static void sdt_eval_grad(void* context, point_t* x, real_t t, vector_t* grad)
{
  // Fetch our Lua table from the registry.
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);

  // Call the grad method.
  lua_getfield(lo->L, -1, "grad");
  lua_pushvalue(lo->L, -2);
  lua_push_point(lo->L, x);
  lua_push_real(lo->L, t);
  lua_call(lo->L, 3, 1);

  if (!lua_is_vector(lo->L, -1))
    luaL_error(lo->L, "grad method did not return a vector.");
  *grad = *lua_to_vector(lo->L, -1);
}

static void sdt_dtor(void* context)
{
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);
  lua_pushnil(lo->L);
  lua_settable(lo->L, LUA_REGISTRYINDEX);
  polymec_free(context);
}

static int sdt_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with a name and methods.");

  // Create a context that houses this object.
  lua_newtable(L);
  int obj_index = num_args + 1;

  sdt_func_vtable vtable = {.dtor = sdt_dtor};

  lua_getfield(L, 1, "name");
  if (!lua_isstring(L, -1))
    return luaL_error(L, "Name must be a string.");
  const char* name = lua_tostring(L, -1);

  lua_getfield(L, 1, "value");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "value must be a method.");
  lua_setfield(L, obj_index, "value");
  vtable.value = sdt_value;

  lua_getfield(L, 1, "grad");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "grad must be a method.");
  lua_setfield(L, obj_index, "grad");
  vtable.eval_grad = sdt_eval_grad;

  // Allocate a context pointer and stuff our object into the registry.
  lua_obj_t* lo = polymec_malloc(sizeof(lua_obj_t));
  lo->L = L;

  // Store the table representing our object in the registry, with 
  // context as a key.
  lua_pushvalue(L, obj_index);
  lua_rawsetp(L, LUA_REGISTRYINDEX, lo);

  // Set up the function.
  sdt_func_t* f = sdt_func_new(name, lo, vtable);
  lua_push_sdt_func(L, f);

  return 1;
}

static int sdt_from_st_funcs(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 3)
    return luaL_error(L, "Arguments must be a name, a signed distance function and its gradient.");
  if (!lua_isstring(L, 1))
    return luaL_error(L, "Argument 1 must be a name.");
  if (!lua_is_st_func(L, 2))
    return luaL_error(L, "Argument 2 must be a time-dependent signed distance function.");
  if (!lua_is_st_func(L, 3))
    return luaL_error(L, "Argument 3 must be the gradient of argument 2.");
  const char* name = lua_tostring(L, 1);
  st_func_t* D = lua_to_st_func(L, 2);
  st_func_t* G = lua_to_st_func(L, 3);
  lua_push_sdt_func(L, sdt_func_from_st_funcs(name, D, G));
  return 1;
}

static lua_module_function sdt_funcs[] = {
  {"new", sdt_new, "sdt_func.new{name = NAME, value = F, grad = gradF} -> new signed distance function."},
  {"from_st_funcs", sdt_from_st_funcs, "sdt_func.from_st_funcs(name, D, G) -> new signed distance function using D for distance and G for gradients."},
  {NULL, NULL, NULL}
};

static int sdt_rename(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a string.");
  sdt_func_rename(f, lua_tostring(L, 2));
  return 0;
}

static int sdt_grad(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument 1 must be a point.");
  point_t* x = lua_to_point(L, 2);
  if (!lua_isnumber(L, 3))
    return luaL_error(L, "Argument 2 must be a time.");
  real_t t = lua_to_real(L, 3);
  vector_t* grad = vector_new(0.0, 0.0, 0.0);
  sdt_func_eval_grad(f, x, t, grad);
  lua_push_vector(L, grad);
  return 1;
}

static int sdt_call(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument 1 must be a point.");
  point_t* x = lua_to_point(L, 2);
  if (!lua_isnumber(L, 3))
    return luaL_error(L, "Argument 2 must be a time.");
  real_t t = lua_to_real(L, 3);
  lua_pushnumber(L, sdt_func_value(f, x, t));
  return 1;
}

static int sdt_project(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument 1 must be a point.");
  if (!lua_isnumber(L, 3))
    return luaL_error(L, "Argument 2 must be a time.");
  point_t* x = lua_to_point(L, 2);
  real_t t = lua_to_real(L, 3);
  point_t* proj_x = point_new(0.0, 0.0, 0.0);
  sdt_func_project(f, x, t, proj_x);
  lua_push_point(L, proj_x);
  return 1;
}

static int sdt_tostring(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  lua_pushfstring(L, "sdt_func '%s'", sdt_func_name(f));
  return 1;
}

static lua_class_method sdt_methods[] = {
  {"rename", sdt_rename, "f:rename(name) -> Renames f to the given name."},
  {"grad", sdt_grad, "f:grad(x) -> Returns the gradient of f at x."},
  {"project", sdt_project, "f:project(x) -> Returns the projection of f to x."},
  {"__call", sdt_call, NULL},
  {"__tostring", sdt_tostring, NULL},
  {NULL, NULL, NULL}
};

static int meshes_uniform_mesh(lua_State* L)
{
  if (!lua_istable(L, 1))
    luaL_error(L, "Argument must be a table with comm, nx, ny, nz, bbox fields.");

  lua_getfield(L, 1, "comm");
  if (!lua_is_mpi_comm(L, -1))
    luaL_error(L, "comm must be an mpi.comm object.");
  MPI_Comm comm = lua_to_mpi_comm(L, -1);

  lua_getfield(L, 1, "nx");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "nx must be a positive integer.");
  int nx = (int)lua_tointeger(L, -1);
  if (nx < 1)
    luaL_error(L, "nx must be positive.");

  lua_getfield(L, 1, "ny");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "ny must be a positive integer.");
  int ny = (int)lua_tointeger(L, -1);
  if (ny < 1)
    luaL_error(L, "ny must be positive.");

  lua_getfield(L, 1, "nz");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "nz must be a positive integer.");
  int nz = (int)lua_tointeger(L, -1);
  if (nz < 1)
    luaL_error(L, "nz must be positive.");

  lua_getfield(L, 1, "bbox");
  if (!lua_is_bbox(L, -1))
    luaL_error(L, "bbox must be a bounding box (bbox).");
  bbox_t* bbox = lua_to_bbox(L, -1);

  int rank = -1;
  lua_getfield(L, 1, "rank");
  if (lua_isinteger(L, -1))
    rank = (int)lua_tointeger(L, -1);

  mesh_t* mesh = NULL;
  if (rank == -1) 
    mesh = create_uniform_mesh(comm, nx, ny, nz, bbox);
  else
    mesh = create_uniform_mesh_on_rank(comm, rank, nx, ny, nz, bbox);

  lua_push_mesh(L, mesh);
  return 1;
}

static real_array_t* get_coordinates(lua_State* L, int index)
{
  real_array_t* array = NULL;
  if (lua_is_array(L, index, LUA_ARRAY_REAL))
    array = lua_to_array(L, index, LUA_ARRAY_REAL);
  else
  {
    array = real_array_new();
    int i = 1;
    while (true)
    {
      lua_rawgeti(L, index, i);
      if (lua_isnil(L, -1))
      {
        lua_pop(L, 1);
        break;
      }
      else
      {
        if (!lua_is_real(L, -1))
          luaL_error(L, "Item %d in table is not a coordinate.");
        real_array_append(array, lua_to_real(L, -1));
        lua_pop(L, 1);
      }
      ++i;
    }
  }
  return array;
}

static int meshes_rectilinear_mesh(lua_State* L)
{
  if (!lua_istable(L, 1))
    luaL_error(L, "Argument must be a table with comm, xs, ys, zs fields.");

  lua_getfield(L, 1, "comm");
  if (!lua_is_mpi_comm(L, -1))
    luaL_error(L, "comm must be an mpi.comm object.");
  MPI_Comm comm = lua_to_mpi_comm(L, -1);

  lua_getfield(L, 1, "xs");
  if (!lua_is_array(L, -1, LUA_ARRAY_REAL) && !lua_istable(L, -1))
    luaL_error(L, "xs must be a table or array of x coordinates.");
  real_array_t* xs = get_coordinates(L, -1);

  lua_getfield(L, 1, "ys");
  if (!lua_is_array(L, -1, LUA_ARRAY_REAL) && !lua_istable(L, -1))
    luaL_error(L, "ys must be a table or array of y coordinates.");
  real_array_t* ys = get_coordinates(L, -1);

  lua_getfield(L, 1, "zs");
  if (!lua_is_array(L, -1, LUA_ARRAY_REAL) && !lua_istable(L, -1))
    luaL_error(L, "zs must be a table or array of z coordinates.");
  real_array_t* zs = get_coordinates(L, -1);

  int rank = -1;
  lua_getfield(L, 1, "rank");
  if (lua_isinteger(L, -1))
    rank = (int)lua_tointeger(L, -1);

  mesh_t* mesh = NULL;
  if (rank == -1) 
  {
    mesh = create_rectilinear_mesh(comm, 
                                   xs->data, (int)xs->size,
                                   ys->data, (int)ys->size,
                                   zs->data, (int)zs->size);
  }
  else
  {
    mesh = create_rectilinear_mesh_on_rank(comm, rank,
                                           xs->data, (int)xs->size,
                                           ys->data, (int)ys->size,
                                           zs->data, (int)zs->size);
  }

  lua_push_mesh(L, mesh);
  return 1;
}

static lua_module_function meshes_funcs[] = {
  {"uniform", meshes_uniform_mesh, "meshes.uniform{comm = COMM, nx = NX, ny = NY, nz = NZ, bbox = BBOX} -> New uniform resolution mesh."},
  {"rectilinear", meshes_rectilinear_mesh, "meshes.rectilinear{comm = COMM, xs = XS, ys = YS, zs = ZS} -> New rectilinear mesh with defined node coordinates."},
  {NULL, NULL, NULL}
};

static lua_module_function points_funcs[] = {
  {NULL, NULL, NULL}
};

//------------------------------------------------------------------------
//                                API 
//------------------------------------------------------------------------

int lua_register_geometry_modules(lua_State* L)
{
  // Core types.
  lua_register_class(L, "sd_func", "A signed distance function.", sd_funcs, sd_methods);
  lua_register_class(L, "sdt_func", "A time-dependent signed distance function.", sdt_funcs, sdt_methods);

  // Register a module of mesh factory methods.
  lua_register_module(L, "meshes", "Functions for generating meshes.", meshes_funcs);

  // Register a module of point factory methods.
  lua_register_module(L, "points", "Functions for generating points.", points_funcs);

  return 0;
}

void lua_push_sd_func(lua_State* L, sd_func_t* f)
{
  lua_push_object(L, "sd_func", f, NULL);
}

bool lua_is_sd_func(lua_State* L, int index)
{
  return lua_is_object(L, index, "sd_func");
}

sd_func_t* lua_to_sd_func(lua_State* L, int index)
{
  return (sd_func_t*)lua_to_object(L, index, "sd_func");
}

void lua_push_sdt_func(lua_State* L, sdt_func_t* f)
{
  lua_push_object(L, "sdt_func", f, NULL);
}

bool lua_is_sdt_func(lua_State* L, int index)
{
  return lua_is_object(L, index, "sdt_func");
}

sdt_func_t* lua_to_sdt_func(lua_State* L, int index)
{
  return (sdt_func_t*)lua_to_object(L, index, "sdt_func");
}

