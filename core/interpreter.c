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

#include "core/interpreter.h"
#include "core/tuple.h"
#include "core/unordered_set.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// This indicates whether the creater/owner of a datum is the Lua interpreter
// or polymec itself.
typedef enum
{
  LUA,
  POLYMEC
} data_origin_t;

// Define a map from variable names to storage.
typedef struct 
{
  void *datum;                 // The datum associated with a variable.
  int size;                    // The size of the datum (if any).
  interpreter_var_type_t type; // The data type.
  void (*dtor)(void*);         // Data destructor.

  // These indicate ownership of data (owner) and metadata (creator).
  data_origin_t owner, creator;
} interpreter_storage_t;

DEFINE_UNORDERED_MAP(interpreter_map, char*, interpreter_storage_t*, string_hash, string_equals)

// This is the destructor for a given storage type within the interpreter.
static int destroy_storage(lua_State* lua)
{
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, -1);
  ASSERT(storage->creator == LUA);
  if ((storage->datum != NULL) && (storage->dtor != NULL) && (storage->owner == LUA))
    (*storage->dtor)(storage->datum);
  storage->datum = NULL;
  return 0;
}

// This destroys variables that have been parsed.
static void destroy_variable(char* key, interpreter_storage_t* value)
{
  free(key);
  if ((value->dtor != NULL) && (value->owner == POLYMEC))
    (*value->dtor)(value->datum);
  value->datum = NULL;
  if (value->creator == POLYMEC)
    free(value);
}

static void destroy_string(char* key)
{
  free(key);
}

static void destroy_table_entry(char* key, void* value)
{
  free(key);
  free(value);
}

static void destroy_table_tuple(char* key, void* value)
{
  free(key);
  real_t* tuple = value;
  real_tuple_free(tuple);
}

// A key-value pair for a metatable.
typedef struct
{
  const char* key;
  lua_CFunction value;
} lua_meta_key_val_t;

// This helper creates a metatable with the given name in the metatable 
// registry, and associates it with the object on top of the Lua stack.
// The table consists of a NULL-terminated array of key-value pairs.
static void set_metatable(lua_State* lua, 
                          const char* table_name,
                          lua_meta_key_val_t table[])
{
  ASSERT(table[0].key != NULL);
  ASSERT(table[0].value != NULL);

  // Register the metatable with the Lua interpreter.
  luaL_newmetatable(lua, table_name);

  int i = 0;
  while (table[i].key != NULL)
  {
    ASSERT(table[i].value != NULL);
    lua_pushcfunction(lua, table[i].value);
    lua_setfield(lua, -2, table[i].key);
    ++i;
  }

  // Set up the destructor for this object.
  lua_pushcfunction(lua, destroy_storage);
  lua_setfield(lua, -2, "__gc");

  // Sets the metatable for the object originally at the top of the stack.
  lua_setmetatable(lua, -2);
}

// This allows us to allocate userdata via the Lua interpreter or 
// after the fact (when lua is NULL).
static interpreter_storage_t* NEW_USER_DATA(lua_State* lua)
{
  interpreter_storage_t* storage;
  if (lua == NULL)
  {
    storage = malloc(sizeof(interpreter_storage_t));
    storage->creator = POLYMEC;
    storage->owner = POLYMEC;
  }
  else
  {
    storage = lua_newuserdata(lua, sizeof(interpreter_storage_t));
    storage->creator = LUA;
    storage->owner = LUA;
  }
  storage->datum = NULL;
  storage->size = 0;
  storage->type = INTERPRETER_TERMINUS;
  storage->dtor = NULL;
  return storage;
}

static interpreter_storage_t* store_string(lua_State* lua, const char* var)
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->datum = string_dup(var);
  storage->type = INTERPRETER_STRING;
  storage->dtor = free;
  return storage;
}

static interpreter_storage_t* store_number(lua_State* lua, real_t var)
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  real_t* dvar = malloc(sizeof(real_t));
  *dvar = var;
  storage->datum = dvar;
  storage->type = INTERPRETER_NUMBER;
  storage->dtor = free;
  return storage;
}

static interpreter_storage_t* store_boolean(lua_State* lua, bool var)
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  bool* bvar = malloc(sizeof(bool));
  *bvar = var;
  storage->datum = bvar;
  storage->type = INTERPRETER_BOOLEAN;
  storage->dtor = free;
  return storage;
}

static interpreter_storage_t* store_point(lua_State* lua, point_t* point)
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->datum = point;
  storage->type = INTERPRETER_POINT;
  storage->dtor = NULL;
  return storage;
}

static interpreter_storage_t* store_pointlist(lua_State* lua, point_t* points, int size)
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->datum = points;
  storage->size = size;
  storage->type = INTERPRETER_POINT_LIST;
  storage->dtor = free;
  return storage;
}

static interpreter_storage_t* store_vector(lua_State* lua, vector_t* vec)
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->datum = vec;
  storage->type = INTERPRETER_VECTOR;
  storage->dtor = NULL;
  return storage;
}

static interpreter_storage_t* store_vectorlist(lua_State* lua, vector_t* vectors, int size)
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->datum = vectors;
  storage->size = size;
  storage->type = INTERPRETER_VECTOR_LIST;
  storage->dtor = free;
  return storage;
}

static interpreter_storage_t* store_boundingbox(lua_State* lua, bbox_t* var)
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->type = INTERPRETER_BOUNDING_BOX;
  storage->datum = var;
  storage->dtor = NULL;
  return storage;
}

static void destroy_mesh(void* mesh)
{
  mesh_free((mesh_t*)mesh);
}

static interpreter_storage_t* store_mesh(lua_State* lua, mesh_t* var)
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->type = INTERPRETER_MESH;
  storage->datum = var;
  storage->dtor = destroy_mesh;
  return storage;
}

static interpreter_storage_t* store_scalar_function(lua_State* lua, st_func_t* var)
{
  ASSERT(st_func_num_comp(var) == 1);
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->type = INTERPRETER_SCALAR_FUNCTION;
  storage->datum = var;
  storage->dtor = NULL;
  return storage;
}

static interpreter_storage_t* store_vector_function(lua_State* lua, st_func_t* var)
{
  ASSERT(st_func_num_comp(var) == 3);
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->type = INTERPRETER_VECTOR_FUNCTION;
  storage->datum = var;
  storage->dtor = NULL;
  return storage;
}

static interpreter_storage_t* store_sym_tensor_function(lua_State* lua, st_func_t* var)
{
  ASSERT(st_func_num_comp(var) == 6);
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->type = INTERPRETER_SYM_TENSOR_FUNCTION;
  storage->datum = var;
  storage->dtor = NULL;
  return storage;
}

static interpreter_storage_t* store_tensor_function(lua_State* lua, st_func_t* var)
{
  ASSERT(st_func_num_comp(var) == 9);
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->type = INTERPRETER_TENSOR_FUNCTION;
  storage->datum = var;
  storage->dtor = NULL;
  return storage;
}

static void destroy_table(void* table)
{
  string_ptr_unordered_map_free((string_ptr_unordered_map_t*)table);
}

static interpreter_storage_t* store_sequence(lua_State* lua, real_t* sequence, int len)
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->datum = sequence;
  storage->type = INTERPRETER_SEQUENCE;
//  storage->dtor = free; // FIXME: Handled by __gc? Why isn't this true for vector/point lists?
  storage->size = len;
  return storage;
}

typedef struct
{
  char** list;
  int len;
} stringlist_t;

static void stringlist_dtor(void* list)
{
  stringlist_t* sl = list;
  for (int i = 0; i < sl->len; ++i)
    free(sl->list[i]);
  free(sl->list);
  free(sl);
}

static interpreter_storage_t* store_stringlist(lua_State* lua, char** list, int len)
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  stringlist_t* sl = malloc(sizeof(stringlist_t));
  sl->list = list;
  sl->len = len;
  storage->datum = sl;
  storage->type = INTERPRETER_STRING_LIST;
  storage->dtor = stringlist_dtor; 
  storage->size = len;
  return storage;
}

static interpreter_storage_t* store_table(lua_State* lua, string_ptr_unordered_map_t* table)
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->datum = table;
  storage->type = INTERPRETER_TABLE;
  storage->dtor = destroy_table;
  return storage;
}

static interpreter_storage_t* store_user_defined(lua_State* lua, void* user_defined, void (*dtor)(void*))
{
  interpreter_storage_t* storage = NEW_USER_DATA(lua);
  storage->type = INTERPRETER_USER_DEFINED;
  storage->datum = user_defined;
  storage->dtor = dtor;
  return storage; 
}

// Interpreter data structure.
struct interpreter_t
{
  // A registry of functions to extend the interpreter.
  int num_functions;
  char* function_names[1024];
  lua_CFunction functions[1024];
  string_array_t* function_docs[1024];

  // A registry of global tables to support OO concepts.
  int num_globals;
  char* global_names[1024];
  string_array_t* global_docs[1024];
  int num_global_methods[1024];
  char* global_method_names[1024][128];
  string_array_t* global_method_docs[1024][128];
  lua_CFunction global_methods[1024][128];

  // The interpreter instance.
  lua_State* lua;

  // The data store.
  interpreter_map_t* store;

  // A list of valid inputs and their types.
  int num_valid_inputs;
  interpreter_validation_t* valid_inputs;

  // A set of pre-existing variables in Lua that will be ignored.
  string_unordered_set_t* preexisting_vars;
};

extern void interpreter_register_default_functions(interpreter_t* interp);

static void interpreter_close_lua(interpreter_t* interp)
{
  if (interp->lua != NULL)
  {
    string_unordered_set_free(interp->preexisting_vars);
    interp->preexisting_vars = NULL;
    lua_close(interp->lua);
    interp->lua = NULL;
  }
}

interpreter_t* interpreter_new(interpreter_validation_t* valid_inputs)
{
  interpreter_t* interp = malloc(sizeof(interpreter_t));
  interp->num_functions = 0;
  interp->num_globals = 0;
  interp->lua = NULL;

  // Add the default functions to the interpreter.
  interpreter_register_default_functions(interp);

  // Initialize the data store.
  interp->store = interpreter_map_new();

  // Copy over the valid inputs.
  if (valid_inputs != NULL)
  {
    int num_valid_inputs = 0;
    while (valid_inputs[num_valid_inputs].type != INTERPRETER_TERMINUS)
      num_valid_inputs++;
    interp->valid_inputs = malloc(num_valid_inputs*sizeof(interpreter_validation_t));
    interp->num_valid_inputs = num_valid_inputs;
  }
  else
  {
    interp->valid_inputs = NULL;
    interp->num_valid_inputs = 0;
  }
  for (int i = 0; i < interp->num_valid_inputs; ++i)
  {
    interp->valid_inputs[i].variable = string_dup(valid_inputs[i].variable);
    interp->valid_inputs[i].type = valid_inputs[i].type;
  }

  interp->preexisting_vars = NULL;

  return interp;
}

void interpreter_free(interpreter_t* interp)
{
  // Put lua away.
  for (int i = 0; i < interp->num_globals; ++i)
  {
    for (int j = 0; j < interp->num_global_methods[i]; ++j)
    {
      free(interp->global_method_names[i][j]);
      if (interp->global_method_docs[i][j] != NULL)
        string_array_free(interp->global_method_docs[i][j]);
    }
    free(interp->global_names[i]);
    if (interp->global_docs[i] != NULL)
      string_array_free(interp->global_docs[i]);
  }
  for (int i = 0; i < interp->num_functions; ++i)
  {
    free(interp->function_names[i]);
    if (interp->function_docs[i] != NULL)
      string_array_free(interp->function_docs[i]);
  }
  for (int i = 0; i < interp->num_valid_inputs; ++i)
    free(interp->valid_inputs[i].variable);
  free(interp->valid_inputs);

  // Make sure we delete our store before closing Lua, since Lua sweeps 
  // all of its variables on exit, so we won't even be able to query the 
  // store to see whether we should be deleting anything.
  interpreter_map_free(interp->store);
  interpreter_close_lua(interp);
  free(interp);
}

void interpreter_register_function(interpreter_t* interp, const char* function_name, int (*function)(lua_State*), string_array_t* doc)
{
  ASSERT(interp->num_functions < 1024);
  interp->function_names[interp->num_functions] = string_dup(function_name);
  interp->functions[interp->num_functions] = function;
  interp->function_docs[interp->num_functions] = doc;
  interp->num_functions++;
}

void interpreter_register_global_table(interpreter_t* interp, const char* table_name, string_array_t* doc)
{
  ASSERT(interp->num_globals < 1024);
  interp->global_names[interp->num_globals] = string_dup(table_name);
  interp->global_docs[interp->num_globals] = doc;
  interp->num_global_methods[interp->num_globals] = 0;
  interp->num_globals++;
}

void interpreter_register_global_method(interpreter_t* interp, const char* table_name, const char* method_name, int (*method)(struct lua_State*), string_array_t* doc)
{
  int global_index = 0;
  while (global_index < interp->num_globals)
  {
    if (!strcmp(table_name, interp->global_names[global_index]))
      break;
    ++global_index;
  }
  if (global_index < interp->num_globals)
  {
    interp->global_method_names[global_index][interp->num_global_methods[global_index]] = string_dup(method_name);
    interp->global_methods[global_index][interp->num_global_methods[global_index]] = method;
    interp->global_method_docs[global_index][interp->num_global_methods[global_index]] = doc;
    interp->num_global_methods[global_index]++;
  }
}

void interpreter_help(interpreter_t* interp, const char* entity, FILE* stream)
{
  ASSERT(entity != NULL);

  // If we were asked to list the available functions, do so.
  if (!strcmp(entity, "list"))
  {
    fprintf(stream, "The following functions are available:\n");
    for (int i = 0; i < interp->num_functions; ++i)
      fprintf(stream, "  %s\n", interp->function_names[i]);
    fprintf(stream, "\nThe following global symbols are available:\n");
    for (int i = 0; i < interp->num_globals; ++i)
      fprintf(stream, "  %s\n", interp->global_names[i]);
    return;
  }

  // Search through registered functions.
  for (int i = 0; i < interp->num_functions; ++i)
  {
    if (!strcmp(entity, interp->function_names[i]))
    {
      // Write the documentation for this function if it exists.
      if (interp->function_docs[i] != NULL)
      {
        for (int j = 0; j < interp->function_docs[i]->size; ++j)
          fprintf(stream, "%s\n", interp->function_docs[i]->data[j]);
      }
      else 
        fprintf(stream, "The function '%s' is currently undocumented.\n", entity);
      return;
    }
  }

  // Search through the global tables/symbols.
  for (int i = 0; i < interp->num_globals; ++i)
  {
    if (!strcmp(entity, interp->global_names[i]))
    {
      // Write the documentation for this global if it exists.
      if (interp->global_docs[i] != NULL)
      {
        for (int j = 0; j < interp->global_docs[i]->size; ++j)
          fprintf(stream, "%s\n", interp->global_docs[i]->data[j]);
      }
      else 
      {
        fprintf(stream, "The global symbol '%s' is currently undocumented.\n", entity);
      }

      // If there are documented symbols are present, write them out.
      for (int j = 0; j < interp->num_global_methods[i]; ++j)
      {
        if (j == 0)
          fprintf(stream, "'%s' contains the following methods:\n", entity);
        if (interp->global_method_docs[i][j] != NULL)
        {
          fprintf(stream, "  %s:\n", interp->global_method_names[i][j]);
          for (int k = 0; k < interp->global_method_docs[i][j]->size; ++k)
            fprintf(stream, "    %s\n", interp->global_method_docs[i][j]->data[k]);
        }
        else
          fprintf(stream, "  %s (undocumented)\n", interp->global_method_names[i][j]);
      }
      return;
    }
  }

  // We didn't find the symbol that was requested.
  fprintf(stream, "No function named '%s' exists.\n", entity);
}

static interpreter_validation_t* interpreter_validation_entry(interpreter_t* interp, const char* key)
{
  for (int i = 0; i < interp->num_valid_inputs; ++i)
  {
    if (!strcmp(key, interp->valid_inputs[i].variable))
      return &interp->valid_inputs[i];
  }
  return NULL;
}

static lua_State* interpreter_open_lua(interpreter_t* interp)
{
  if (interp->lua != NULL)
    interpreter_close_lua(interp);

  // Initialize the Lua interpreter.
  interp->lua = luaL_newstate();
  ASSERT(interp->lua != NULL);
  luaL_openlibs(interp->lua);

  // Add whatever other functions we've had registered.
  for (int i = 0; i < interp->num_functions; ++i)
    lua_register(interp->lua, (const char*)interp->function_names[i], interp->functions[i]);

  // Add whatever global variables we've had registered (and whatever
  // methods we've attached to them).
  for (int i = 0; i < interp->num_globals; ++i)
  {
    lua_createtable(interp->lua, 16, 16);
    lua_setglobal(interp->lua, interp->global_names[i]);
    lua_getglobal(interp->lua, interp->global_names[i]);
    int table_index = -2;
    for (int j = 0; j < interp->num_global_methods[i]; ++j)
    {
      lua_pushcfunction(interp->lua, interp->global_methods[i][j]);
      lua_setfield(interp->lua, table_index, interp->global_method_names[i][j]);
    }
  }

  // Before we go adding our own variables to the interpreter, make a 
  // list of what's there so we can ignore it.
  // FIXME: This is a hack, but it's easier than learning about 
  // FIXME: lua's notions of environments, etc.
  ASSERT(interp->preexisting_vars == NULL);
  interp->preexisting_vars = string_unordered_set_new();
  lua_pushglobaltable(interp->lua);
  ASSERT(lua_istable(interp->lua, -1));
  lua_pushnil(interp->lua); 
  while (lua_next(interp->lua, -2)) // Traverse globals table.
  {
    const char* key = lua_tostring(interp->lua, -2);
    string_unordered_set_insert(interp->preexisting_vars, (char*)key);
    lua_pop(interp->lua, 1);
  }
  lua_pop(interp->lua, 1);

  return interp->lua;
}

static void interpreter_store_chunk_contents(interpreter_t* interp)
{
  ASSERT(interp->lua != NULL);
  lua_State* lua = interp->lua;

  // Get the global variables table from the registry.
  lua_pushglobaltable(lua);
//  lua_rawgeti(lua, LUA_REGISTRYINDEX, LUA_RIDX_GLOBALS);
  ASSERT(lua_istable(lua, -1));

  // Traverse the table.
  lua_pushnil(lua); 
  while (lua_next(lua, -2)) // Traverse globals table.
  {
    // key is at index -2, value is at index -1.
    static const int key_index = -2;
    static const int val_index = -1;
    const char* key = lua_tostring(lua, key_index);

    // Check to see if this variable existed in Lua before we 
    // parsed our inputs.
    bool preexisting_var = string_unordered_set_contains(interp->preexisting_vars, (char*)key);
    bool skip_this_var = false;

    // Does the key appear in our validation table?
    interpreter_validation_t* entry = interpreter_validation_entry(interp, key);
    if (entry != NULL)
    {
      // We must validate this variable against its allowed type.
      if ((entry->type == INTERPRETER_STRING) && !lua_isstring(lua, val_index))
      {
        if (preexisting_var)
          skip_this_var = true;
        else
          polymec_error("Type error: %s must be a string.", key);
      }
      else if ((entry->type == INTERPRETER_NUMBER) && !lua_isnumber(lua, val_index))
      {
        if (preexisting_var)
          skip_this_var = true;
        else
          polymec_error("Type error: %s must be a number.", key);
      }
      else if (entry->type == INTERPRETER_POINT)
      {
        if (!lua_ispoint(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a point.", key);
        }
      }
      else if (entry->type == INTERPRETER_POINT_LIST)
      {
        if (!lua_ispointlist(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a list of points.", key);
        }
      }
      else if (entry->type == INTERPRETER_VECTOR)
      {
        if (!lua_isvector(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a vector.", key);
        }
      }
      else if (entry->type == INTERPRETER_VECTOR_LIST)
      {
        if (!lua_isvectorlist(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a list of vectors.", key);
        }
      }
      else if (entry->type == INTERPRETER_BOUNDING_BOX)
      {
        if (!lua_isboundingbox(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a bounding box.", key);
        }
      }
      else if (entry->type == INTERPRETER_MESH)
      {
        if (!lua_isuserdata(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a mesh.", key);
        }
        interpreter_storage_t* var = (void*)lua_topointer(lua, val_index);
        if (var->type != INTERPRETER_MESH)
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a mesh.", key);
        }
      }
      else if (entry->type == INTERPRETER_SCALAR_FUNCTION)
      {
        if (!lua_isuserdata(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a scalar-valued function.", key);
        }
        interpreter_storage_t* var = (void*)lua_topointer(lua, val_index);
        if (var->type != INTERPRETER_SCALAR_FUNCTION)
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a scalar-valued function.", key);
        }
      }
      else if (entry->type == INTERPRETER_VECTOR_FUNCTION)
      {
        if (!lua_isuserdata(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a vector-valued function.", key);
        }
        interpreter_storage_t* var = (void*)lua_topointer(lua, val_index);
        if (var->type != INTERPRETER_VECTOR_FUNCTION)
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a vector-valued function.", key);
        }
      }
      else if (entry->type == INTERPRETER_SYM_TENSOR_FUNCTION)
      {
        if (!lua_isnumber(lua, val_index))
        {
          if (!lua_isuserdata(lua, val_index))
          {
            if (preexisting_var)
              skip_this_var = true;
            else
              polymec_error("Type error: %s must be a symmetric-tensor-valued function.", key);
          }
          interpreter_storage_t* var = (void*)lua_topointer(lua, val_index);
          if (var->type != INTERPRETER_SYM_TENSOR_FUNCTION)
          {
            if (preexisting_var)
              skip_this_var = true;
            else
              polymec_error("Type error: %s must be a symmetric-tensor-valued function.", key);
          }
        }
      }
      else if (entry->type == INTERPRETER_TENSOR_FUNCTION)
      {
        if (!lua_isnumber(lua, val_index))
        {
          if (!lua_isuserdata(lua, val_index))
          {
            if (preexisting_var)
              skip_this_var = true;
            else
              polymec_error("Type error: %s must be a tensor-valued function.", key);
          }
          interpreter_storage_t* var = (void*)lua_topointer(lua, val_index);
          if (var->type != INTERPRETER_TENSOR_FUNCTION)
          {
            if (preexisting_var)
              skip_this_var = true;
            else
              polymec_error("Type error: %s must be a tensor-valued function.", key);
          }
        }
      }
      else if ((entry->type == INTERPRETER_SEQUENCE) && !lua_issequence(lua, val_index))
      {
        if (preexisting_var)
          skip_this_var = true;
        else
          polymec_error("Type error: %s must be a sequence of numbers.", key);
      }
      else if ((entry->type == INTERPRETER_STRING_LIST) && !lua_isstringlist(lua, val_index))
      {
        if (preexisting_var)
          skip_this_var = true;
        else
          polymec_error("Type error: %s must be a list of strings.", key);
      }
      else if ((entry->type == INTERPRETER_TABLE) && !lua_istable(lua, val_index))
      {
        if (preexisting_var)
          skip_this_var = true;
        else
          polymec_error("Type error: %s must be a table mapping strings to objects.", key);
      }

      // Skip empty tables, as they cause problems.
      else if (lua_istable(lua, val_index))
      {
        // Count the entries. Neither lua_rawlen() nor luaL_len() 
        // works for non-sequence tables.
        int tab_index = lua_absindex(lua, val_index);
        int num_entries = 0;
        lua_pushnil(lua);
        while (lua_next(lua, tab_index))
        {
          ++num_entries;
          lua_pop(lua, 1);
        }
        if (num_entries == 0)
          skip_this_var = true;
      }
    }

    // We skip C functions altogether.
    if (lua_iscfunction(lua, val_index))
      skip_this_var = true;

    // Skip this variable if we need to.
    if (skip_this_var)
    {
      lua_pop(lua, 1);
      continue;
    }

    // Now store the present lua variable in appropriate C storage.
    interpreter_storage_t* var = NULL;
    if (lua_isnumber(lua, val_index))
      var = store_number(NULL, (real_t)lua_tonumber(lua, val_index));
    else if (lua_isboolean(lua, val_index))
      var = store_boolean(NULL, lua_toboolean(lua, val_index));
    else if (lua_isstring(lua, val_index))
      var = store_string(NULL, lua_tostring(lua, val_index));
    else if (lua_issequence(lua, val_index))
    {
      // Sequences can be interpreted in many ways. Are we asked to 
      // interpret this sequence as something other than just a sequence?
      int len;
      real_t* seq = lua_tosequence(lua, val_index, &len);
      if ((entry == NULL) || (entry->type == INTERPRETER_SEQUENCE))
      {
        var = store_sequence(NULL, seq, len);
      }
      else if (entry->type == INTERPRETER_POINT)
      {
        point_t* p = point_new(seq[0], seq[1], seq[2]);
        var = store_point(NULL, p);
      }
      else if (entry->type == INTERPRETER_VECTOR)
      {
        vector_t* v = vector_new(seq[0], seq[1], seq[2]);
        var = store_vector(NULL, v);
      }

      // Clean up if necessary.
      if ((entry != NULL) && (entry->type != INTERPRETER_SEQUENCE))
        free(seq);
    }
    else if (lua_isstringlist(lua, val_index))
    {
      int len;
      char** list = lua_tostringlist(lua, val_index, &len);
      var = store_stringlist(NULL, list, len);
    }
    else if (lua_isboundingbox(lua, val_index))
    {
      var = (void*)lua_topointer(lua, val_index);
    }
    else if (lua_ispointlist(lua, val_index))
    {
      // Point lists and vector lists are interchangeable, so we will 
      // simply assume that either is a point list.
      if (lua_istable(lua, val_index))
      {
        int size;
        point_t* pointlist = lua_topointlist(lua, val_index, &size);
        var = store_pointlist(NULL, pointlist, size);
      }
      else
        var = (void*)lua_topointer(lua, val_index);
    }
    else if (lua_istable(lua, val_index))
    {
      // Before we do anything, we validate the table.
      interpreter_var_type_t value_type = INTERPRETER_TERMINUS;
      static const char* type_names[] = {"string", "number", "boolean", "sequence", "mesh", "function"};
      // Traverse this table and make sure its values are all of one type.
      lua_pushnil(lua);
      while (lua_next(lua, -2))
      {
        // Key is at index -2, value is at -1.
        static const int key_index = -2;
        static const int val_index = -1;
        if (!lua_isstring(lua, key_index) && !lua_isnumber(lua, key_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a list of objects or a table mapping strings to objects.", key);
        }
        if (!lua_isnumber(lua, val_index) && 
            !lua_isboolean(lua, val_index) && 
            !lua_isstring(lua, val_index) && 
            !lua_issequence(lua, val_index) && 
            !lua_isuserdata(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a list of objects or a table mapping strings to objects.", key);
        }
        if (skip_this_var)
        {
          lua_pop(lua, 1);
          continue;
        }

        // Convert our key to a string. If it's an integer, it will be 
        // converted to its string representation.
        char* tkey = (char*)lua_tostring(lua, key_index);
        void* tval = lua_touserdata(lua, val_index);
        if (tval == NULL)
        {
          // We don't need to do any further validation on non-pointer 
          // objects.
          lua_pop(lua, 1);
          continue;
        }

        // No tables of tables allowed!
        if (preexisting_var)
          skip_this_var = true;
        else
        {
          interpreter_storage_t* tvar = tval;
          if (tvar->type == INTERPRETER_TABLE)
            polymec_error("Value error: Key '%s' in table %s stores a table (not allowed!)", tkey, key); 
          if (value_type == INTERPRETER_TERMINUS)
            value_type = tvar->type;
          else if (tvar->type != value_type)
          {
            polymec_error("Value error: Key '%s' in table %s stores a %s (should be %s)", 
                tkey, key, type_names[tvar->type], type_names[value_type]);
          }
        }

        // Removes value from stack.
        lua_pop(lua, 1);
      }

      // If we need to skip the table/sequence, do so.
      if (skip_this_var)
      {
        lua_pop(lua, 1);
        continue;
      }

      // Now traverse the table again and create a C data structure.
      string_ptr_unordered_map_t* table = string_ptr_unordered_map_new();
      lua_pushnil(lua);
      while (lua_next(lua, -2))
      {
        // Key is at index -2, value is at -1.
        static const int key_index = -2;
        static const int val_index = -1;
        char* tkey = (char*)lua_tostring(lua, key_index);
        if (lua_isnumber(lua, val_index))
        {
          real_t* var = malloc(sizeof(real_t));
          *var = (real_t)lua_tonumber(lua, val_index);
          string_ptr_unordered_map_insert_with_kv_dtor(table, tkey, var, destroy_table_entry);
        }
        else if (lua_isboolean(lua, val_index))
        {
          bool* var = malloc(sizeof(bool));
          *var = lua_toboolean(lua, val_index);
          string_ptr_unordered_map_insert_with_kv_dtor(table, tkey, var, destroy_table_entry);
        }
        else if (lua_isstring(lua, val_index))
        {
          const char* var = lua_tostring(lua, val_index);
          string_ptr_unordered_map_insert_with_kv_dtor(table, tkey, string_dup(var), destroy_table_entry);
        }
        else if (lua_issequence(lua, val_index))
        {
          // We store sequences as tuples.
          int len;
          real_t* var = lua_tosequence(lua, val_index, &len);
          real_t* tuple = real_tuple_new(len);
          memcpy(tuple, var, sizeof(real_t) * len);
          free(var);
          string_ptr_unordered_map_insert_with_kv_dtor(table, tkey, tuple, destroy_table_tuple);
        }
        else if (lua_isuserdata(lua, val_index))
        {
          void* tval = (void*)lua_topointer(lua, val_index);
          interpreter_storage_t* tvar = (interpreter_storage_t*)tval;
          string_ptr_unordered_map_insert_with_kv_dtor(table, tkey, tvar->datum, destroy_table_entry);
        }

        // Removes value from stack.
        lua_pop(lua, 1);
      }
      var = store_table(NULL, table);
    }
    else 
    {
      // We're out of types, so we hope this one's a user data.
      if (!lua_isuserdata(lua, val_index))
      {
        // Skip this variable if it existed before and is the wrong type.
        if (preexisting_var)
        {
          lua_pop(lua, 1);
          continue;
        }
      }
      var = (void*)lua_topointer(lua, val_index);
    }

    interpreter_map_insert_with_k_dtor(interp->store, string_dup(key), var, destroy_string); 

    // Removes value from stack -- key is kept for next iteration.
    lua_pop(lua, 1);
  }
  lua_pop(lua, 1); // Pops the globals table.
}

void interpreter_parse_string(interpreter_t* interp, char* input_string)
{
  // Clear the current data store.
  interpreter_map_clear(interp->store);

  // (Re-)initialize the Lua interpreter.
  interpreter_open_lua(interp);

  // Parse the input and execute it.
  int errors = luaL_dostring(interp->lua, input_string);
  if (errors != 0)
    polymec_error(lua_tostring(interp->lua, -1));

  // Move the data from the chunk to the data store.
  interpreter_store_chunk_contents(interp);
}

void interpreter_parse_file(interpreter_t* interp, char* input_file)
{
  log_detail("interpreter: Looking for input in file '%s'...", input_file);

  // Clear the current data store.
  interpreter_map_clear(interp->store);

  // (Re-)initialize the Lua interpreter.
  interpreter_open_lua(interp);

  // Load the file and execute it.
  int errors = luaL_dofile(interp->lua, input_file);
  if (errors != 0)
    polymec_error(lua_tostring(interp->lua, -1));

  // Move the data from the chunk to the data store.
  interpreter_store_chunk_contents(interp);
}

bool interpreter_contains(interpreter_t* interp, const char* variable, interpreter_var_type_t type)
{
  ASSERT(type != INTERPRETER_TERMINUS);
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)variable);
  if (storage == NULL)
    return false;
  return ((*storage)->type == type);
}

char* interpreter_get_string(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_STRING)
    return NULL;
  (*storage)->owner = POLYMEC;
  return (char*)((*storage)->datum);
}

void interpreter_set_string(interpreter_t* interp, const char* name, const char* value)
{
  interpreter_storage_t* storage = store_string(NULL, value);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

real_t interpreter_get_number(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return -FLT_MAX;
  if ((*storage)->type != INTERPRETER_NUMBER)
    return -FLT_MAX;
  return *((real_t*)(*storage)->datum);
}

void interpreter_set_number(interpreter_t* interp, const char* name, real_t value)
{
  interpreter_storage_t* storage = store_number(NULL, value);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

bool interpreter_get_boolean(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return false;
  if ((*storage)->type != INTERPRETER_BOOLEAN)
    return false;
  return *((bool*)(*storage)->datum);
}

void interpreter_set_boolean(interpreter_t* interp, const char* name, bool value)
{
  interpreter_storage_t* storage = store_boolean(NULL, value);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

point_t* interpreter_get_point(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_POINT)
    return NULL;
  (*storage)->owner = POLYMEC;
  return (point_t*)((*storage)->datum);
}

void interpreter_set_point(interpreter_t* interp, const char* name, point_t* value)
{
  interpreter_storage_t* storage = store_point(NULL, value);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

point_t* interpreter_get_pointlist(interpreter_t* interp, const char* name, int* num_points)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_POINT_LIST)
    return NULL;
  (*storage)->owner = POLYMEC;
  *num_points = (*storage)->size;
  return (point_t*)((*storage)->datum);
}

void interpreter_set_pointlist(interpreter_t* interp, const char* name, point_t* points, int num_points)
{
  interpreter_storage_t* storage = store_pointlist(NULL, points, num_points);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

vector_t* interpreter_get_vector(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_VECTOR)
    return NULL;
  (*storage)->owner = POLYMEC;
  return (vector_t*)((*storage)->datum);
}

void interpreter_set_vector(interpreter_t* interp, const char* name, vector_t* value)
{
  interpreter_storage_t* storage = store_vector(NULL, value);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

vector_t* interpreter_get_vectorlist(interpreter_t* interp, const char* name, int* num_vectors)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_VECTOR_LIST)
    return NULL;
  (*storage)->owner = POLYMEC;
  *num_vectors = (*storage)->size;
  return (vector_t*)((*storage)->datum);
}

void interpreter_set_vectorlist(interpreter_t* interp, const char* name, vector_t* vectors, int num_vectors)
{
  interpreter_storage_t* storage = store_vectorlist(NULL, vectors, num_vectors);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

bbox_t* interpreter_get_boundingbox(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_BOUNDING_BOX)
    return NULL;
  (*storage)->owner = POLYMEC;
  return (bbox_t*)((*storage)->datum);
}

void interpreter_set_boundingbox(interpreter_t* interp, const char* name, bbox_t* value)
{
  interpreter_storage_t* storage = store_boundingbox(NULL, value);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

mesh_t* interpreter_get_mesh(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_MESH)
    return NULL;
  (*storage)->owner = POLYMEC;
  return (mesh_t*)((*storage)->datum);
}

void interpreter_set_mesh(interpreter_t* interp, const char* name, mesh_t* value)
{
  interpreter_storage_t* storage = store_mesh(NULL, value);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

st_func_t* interpreter_get_scalar_function(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_SCALAR_FUNCTION)
    return NULL;
  st_func_t* func = (st_func_t*)((*storage)->datum);
  ASSERT(st_func_num_comp(func) == 1);
  (*storage)->owner = POLYMEC;
  return func;
}

void interpreter_set_scalar_function(interpreter_t* interp, const char* name, st_func_t* value)
{
  interpreter_storage_t* storage = store_scalar_function(NULL, value);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

st_func_t* interpreter_get_vector_function(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_VECTOR_FUNCTION)
    return NULL;
  st_func_t* func = (st_func_t*)((*storage)->datum);
  ASSERT(st_func_num_comp(func) == 3);
  (*storage)->owner = POLYMEC;
  return func;
}

void interpreter_set_vector_function(interpreter_t* interp, const char* name, st_func_t* value)
{
  interpreter_storage_t* storage = store_vector_function(NULL, value);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

st_func_t* interpreter_get_sym_tensor_function(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type == INTERPRETER_NUMBER)
  {
    real_t F = *(real_t*)((*storage)->datum);
    real_t v[6] = {F, 0.0, 0.0, F, 0.0, F};
    return constant_st_func_new(6, v);
  }
  if ((*storage)->type != INTERPRETER_SYM_TENSOR_FUNCTION)
    return NULL;
  st_func_t* func = (st_func_t*)((*storage)->datum);
  ASSERT(st_func_num_comp(func) == 6);
  (*storage)->owner = POLYMEC;
  return func;
}

void interpreter_set_sym_tensor_function(interpreter_t* interp, const char* name, st_func_t* value)
{
  interpreter_storage_t* storage = store_sym_tensor_function(NULL, value);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

st_func_t* interpreter_get_tensor_function(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type == INTERPRETER_NUMBER)
  {
    real_t F = *(real_t*)((*storage)->datum);
    real_t v[9] = {F, 0.0, 0.0, 0.0, F, 0.0, 0.0, 0.0, F};
    return constant_st_func_new(9, v);
  }
  if ((*storage)->type != INTERPRETER_TENSOR_FUNCTION)
    return NULL;
  st_func_t* func = (st_func_t*)((*storage)->datum);
  ASSERT(st_func_num_comp(func) == 9);
  (*storage)->owner = POLYMEC;
  return func;
}

void interpreter_set_tensor_function(interpreter_t* interp, const char* name, st_func_t* value)
{
  interpreter_storage_t* storage = store_tensor_function(NULL, value);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

string_ptr_unordered_map_t* interpreter_get_table(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_TABLE)
    return NULL;
  (*storage)->owner = POLYMEC;
  return (string_ptr_unordered_map_t*)((*storage)->datum);
}

void interpreter_set_table(interpreter_t* interp, const char* name, string_ptr_unordered_map_t* value)
{
  interpreter_storage_t* storage = store_table(NULL, value);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

real_t* interpreter_get_sequence(interpreter_t* interp, const char* name, int* len)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_SEQUENCE)
    return NULL;
  (*storage)->owner = POLYMEC;
  *len = (*storage)->size;
  return (real_t*)((*storage)->datum);
}

void interpreter_set_sequence(interpreter_t* interp, const char* name, real_t* sequence, int len)
{
  interpreter_storage_t* storage = store_sequence(NULL, sequence, len);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

char** interpreter_get_stringlist(interpreter_t* interp, const char* name, int* len)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_STRING_LIST)
    return NULL;
  (*storage)->owner = POLYMEC;
  stringlist_t* sl = (*storage)->datum;
  ASSERT((*storage)->size == sl->len);
  *len = sl->len;
  return sl->list;
}

void interpreter_set_stringlist(interpreter_t* interp, const char* name, char** list, int len)
{
  interpreter_storage_t* storage = store_stringlist(NULL, list, len);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

void* interpreter_get_user_defined(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_USER_DEFINED)
    return NULL;
  (*storage)->owner = POLYMEC;
  return (*storage)->datum;
}

void interpreter_set_user_defined(interpreter_t* interp, const char* name, void* value, void (*dtor)(void*))
{
  interpreter_storage_t* storage = store_user_defined(NULL, value, dtor);
  interpreter_map_insert_with_kv_dtor(interp->store, string_dup(name), storage, destroy_variable);
}

//------------------------------------------------------------------------
//                          Lua helpers.
//------------------------------------------------------------------------

bool lua_issequence(struct lua_State* lua, int index)
{
  index = lua_absindex(lua, index);
  if (lua_istable(lua, index))
  {
    lua_pushnil(lua);
    while (lua_next(lua, index))
    {
      // Key is at index -2, value is at -1.
      static const int key_index = -2;
      static const int val_index = -1;
      bool key_is_number = lua_isnumber(lua, key_index);
      bool val_is_number = lua_isnumber(lua, val_index);
      lua_pop(lua, 1);
      if (!key_is_number || !val_is_number)
      {
        lua_pop(lua, 1);
        return false;
      }
    }
    return true;
  }
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_SEQUENCE);
}

real_t* lua_tosequence(struct lua_State* lua, int index, int* len)
{
  if (!lua_issequence(lua, index))
    return NULL;
  index = lua_absindex(lua, index);
  if (lua_istable(lua, index))
  {
    *len = lua_rawlen(lua, index);
    real_t* seq = malloc(sizeof(real_t)*(*len));
    for (int i = 1; i <= *len; ++i)
    {
      lua_pushinteger(lua, (lua_Integer)i);
      lua_gettable(lua, index);
      seq[i-1] = (real_t)lua_tonumber(lua, -1);
      lua_pop(lua, 1);
    }
    return seq;
  }
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_SEQUENCE)
  {
    *len = storage->size;
    return (real_t*)storage->datum;
  }
  else
    return NULL;
}

static int sequence_tostring(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_SEQUENCE);
  real_t* data = var->datum;
  char* str = malloc(sizeof(char) * 18 * var->size);
  str[0] = '{';
  int offset = 1;
  for (int i = 0; i < var->size; ++i)
  {
    char stri[19];
    if (i < (var->size-1))
      sprintf(stri, "%g, ", data[i]);
    else
      sprintf(stri, "%g}", data[i]);
    strcpy(&str[offset], stri);
    offset += strlen(stri);
  }
  lua_pushstring(lua, str);
  free(str);
  return 1;
}

static int sequence_concat(lua_State* lua)
{
  interpreter_storage_t* var1 = (void*)lua_topointer(lua, -1);
  ASSERT(var1->type == INTERPRETER_SEQUENCE);
  real_t* data1 = var1->datum;
  interpreter_storage_t* var2 = (void*)lua_topointer(lua, -2);
  ASSERT(var2->type == INTERPRETER_SEQUENCE);
  real_t* data2 = var2->datum;

  int len = var1->size + var2->size;
  real_t* concat_data = malloc(sizeof(real_t) * len);
  memcpy(concat_data, data1, sizeof(real_t) * var1->size);
  memcpy(&concat_data[var1->size], data2, sizeof(real_t) * var2->size);
  lua_pushsequence(lua, concat_data, len);
  return 1;
}

static int sequence_len(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_SEQUENCE);
  lua_pushinteger(lua, var->size);
  return 1;
}

void lua_pushsequence(struct lua_State* lua, real_t* sequence, int len)
{
  // Bundle it up and store it in the given variable.
  store_sequence(lua, sequence, len);
  lua_meta_key_val_t metatable[] = 
    {{"__tostring", sequence_tostring},
     {"__concat", sequence_concat},
     {"__len", sequence_len},
     {NULL, NULL}};
  set_metatable(lua, "sequence_metatable", metatable);
}

bool lua_isstringlist(struct lua_State* lua, int index)
{
  index = lua_absindex(lua, index);
  if (lua_istable(lua, index))
  {
    lua_pushnil(lua);
    while (lua_next(lua, index))
    {
      // Key is at index -2, value is at -1.
      static const int key_index = -2;
      static const int val_index = -1;
      bool key_is_number = lua_isnumber(lua, key_index);
      bool val_is_string = lua_isstring(lua, val_index);
      lua_pop(lua, 1);
      if (!key_is_number || !val_is_string)
      {
        lua_pop(lua, 1);
        return false;
      }
    }
    return true;
  }
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_STRING_LIST);
}

char** lua_tostringlist(struct lua_State* lua, int index, int* len)
{
  if (!lua_isstringlist(lua, index))
    return NULL;
  index = lua_absindex(lua, index);
  if (lua_istable(lua, index))
  {
    *len = lua_rawlen(lua, index);
    char** list = malloc(sizeof(char*)*(*len));
    for (int i = 1; i <= *len; ++i)
    {
      lua_pushinteger(lua, (lua_Integer)i);
      lua_gettable(lua, index);
      list[i-1] = string_dup(lua_tostring(lua, -1));
      lua_pop(lua, 1);
    }
    return list;
  }
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_STRING_LIST)
  {
    stringlist_t* sl = storage->datum;
    ASSERT(sl->len == *len);
    *len = sl->len;
    return sl->list;
  }
  else
    return NULL;
}

static int stringlist_tostring(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_SEQUENCE);
  stringlist_t* sl = var->datum;
  char** data = sl->list;
  int repr_len = 0;
  for (int i = 0; i < sl->len; ++i)
    repr_len += strlen(data[i]) + 10;
  char* str = malloc(sizeof(char) * repr_len);
  str[0] = '{';
  int offset = 1;
  for (int i = 0; i < var->size; ++i)
  {
    char stri[strlen(data[i])+10];
    if (i < (var->size-1))
      sprintf(stri, "\'%s\', ", data[i]);
    else
      sprintf(stri, "'%s'}", data[i]);
    strcpy(&str[offset], stri);
    offset += strlen(stri);
  }
  lua_pushstring(lua, str);
  free(str);
  return 1;
}

static int stringlist_concat(lua_State* lua)
{
  interpreter_storage_t* var1 = (void*)lua_topointer(lua, -1);
  ASSERT(var1->type == INTERPRETER_STRING_LIST);
  char** data1 = var1->datum;
  interpreter_storage_t* var2 = (void*)lua_topointer(lua, -2);
  ASSERT(var2->type == INTERPRETER_STRING_LIST);
  char** data2 = var2->datum;

  int len = var1->size + var2->size;
  char** concat_data = malloc(sizeof(char*) * len);
  for (int i = 0; i < var1->size; ++i)
    concat_data[i] = string_dup(data1[i]);
  for (int i = 0; i < var2->size; ++i)
    concat_data[i+var1->size] = string_dup(data2[i]);
  lua_pushstringlist(lua, concat_data, len);
  return 1;
}

static int stringlist_len(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_STRING_LIST);
  lua_pushinteger(lua, var->size);
  return 1;
}

void lua_pushstringlist(struct lua_State* lua, char** list, int len)
{
  // Bundle it up and store it in the given variable.
  store_stringlist(lua, list, len);
  lua_meta_key_val_t metatable[] = 
    {{"__tostring", stringlist_tostring},
     {"__concat", stringlist_concat},
     {"__len", stringlist_len},
     {NULL, NULL}};
  set_metatable(lua, "stringlist_metatable", metatable);
}

bool lua_ispoint(struct lua_State* lua, int index)
{
  // A sequence with 3 numbers in it will work.
  if (lua_issequence(lua, index) && (lua_rawlen(lua, index) == 3))
    return true;
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_POINT);
}

point_t* lua_topoint(struct lua_State* lua, int index)
{
  if (!lua_ispoint(lua, index))
    return NULL;
  if (lua_istable(lua, index))
  {
    int len;
    real_t* seq = lua_tosequence(lua, index, &len);
    point_t* p = point_new(seq[0], seq[1], seq[2]);
    free(seq);
    return p;
  }
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_POINT)
    return (point_t*)storage->datum;
  else if (storage->type == INTERPRETER_SEQUENCE)
  {
    real_t* seq = storage->datum;
    return point_new(seq[0], seq[1], seq[2]);
  }
  else
    return NULL;
}

static int point_tostring(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_POINT);
  point_t* data = var->datum;
  char str[80];
  sprintf(str, "{%g, %g, %g}", data->x, data->y, data->z);
  lua_pushstring(lua, str);
  return 1;
}

static int point_index(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_POINT);
  if (!lua_isnumber(lua, -1))
    return luaL_error(lua, "Non-numeric index given for point.");
  int index = (int)lua_tonumber(lua, -1);
  if ((index < 1) || (index > 3))
    return luaL_error(lua, "Invalid index for point: %d", index);
  point_t* data = var->datum;
  real_t comp = (index == 1) ? data->x : (index == 2) ? data->y : data->z;
  lua_pushnumber(lua, comp);
  return 1;
}

static int point_mul(lua_State* lua)
{
  interpreter_storage_t* var;
  int factor;
  if (lua_isnumber(lua, -1))
  {
    factor = (int)lua_tonumber(lua, -1);
    var = (void*)lua_topointer(lua, -2);
  }
  else
  {
    factor = (int)lua_tonumber(lua, -2);
    var = (void*)lua_topointer(lua, -1);
  }
  ASSERT(var->type == INTERPRETER_POINT);

  point_t* data = var->datum;
  point_t* ptlist = malloc(sizeof(point_t) * factor);
  for (int i = 0; i < factor; ++i)
    ptlist[i] = *data;
  lua_pushpointlist(lua, ptlist, factor);
  return 1;
}

void lua_pushpoint(struct lua_State* lua, point_t* point)
{
  // Bundle it up and store it in the given variable.
  store_point(lua, point);
  lua_meta_key_val_t metatable[] = 
    {{"__tostring", point_tostring},
     {"__index", point_index},
     {"__mul", point_mul},
     {NULL, NULL}};
  set_metatable(lua, "point_metatable", metatable);
}

bool lua_ispointlist(struct lua_State* lua, int index)
{
  index = lua_absindex(lua, index);
  if (lua_istable(lua, index)) // A table with points in it will work.
  {
    size_t len = lua_rawlen(lua, index);
    if (len == 0) 
      return false;
    for (size_t i = 1; i <= len; ++i)
    {
      lua_pushinteger(lua, (lua_Integer)i);
      lua_gettable(lua, index);
      bool is_point = lua_ispoint(lua, -1);
      lua_pop(lua, 1);
      if (!is_point) 
        return false;
    }
    return true;
  }
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_POINT_LIST);
}

point_t* lua_topointlist(struct lua_State* lua, int index, int* size)
{
  if (!lua_ispointlist(lua, index))
    return NULL;
  index = lua_absindex(lua, index);
  if (lua_istable(lua, index))
  {
    *size = (int)lua_rawlen(lua, index);
    point_t* points = malloc(sizeof(point_t) * (*size));
    for (int i = 0; i < *size; ++i)
    {
      lua_pushinteger(lua, (lua_Integer)(i+1));
      lua_gettable(lua, index);
      point_t* p = lua_topoint(lua, -1);
      lua_pop(lua, 1);
      point_copy(&points[i], p);
    }
    return points;
  }
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_POINT_LIST)
  {
    *size = storage->size;
    return (point_t*)storage->datum;
  }
  else
    return NULL;
}

static int pointlist_tostring(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_POINT_LIST);
  point_t* data = var->datum;
  char* str = malloc(sizeof(char) * 80 * var->size);
  str[0] = '{';
  int offset = 1;
  for (int i = 0; i < var->size; ++i)
  {
    char stri[81];
    if (i < (var->size-1))
      sprintf(stri, "{%g, %g, %g}, ", data[i].x, data[i].y, data[i].z);
    else
      sprintf(stri, "{%g, %g, %g}}", data[i].x, data[i].y, data[i].z);
    strcpy(&str[offset], stri);
    offset += strlen(stri);
  }
  lua_pushstring(lua, str);
  free(str);
  return 1;
}

static int pointlist_concat(lua_State* lua)
{
  interpreter_storage_t* var1 = (void*)lua_topointer(lua, -1);
  ASSERT(var1->type == INTERPRETER_POINT_LIST);
  point_t* data1 = var1->datum;
  if (!lua_ispointlist(lua, -2))
    polymec_error("Only lists of points can be concatenated to point lists.");
  int num_points;
  point_t* points = lua_topointlist(lua, -2, &num_points);

  int len = var1->size + num_points;
  point_t* concat_data = malloc(sizeof(point_t) * len);
  memcpy(concat_data, data1, sizeof(point_t) * var1->size);
  memcpy(&concat_data[var1->size], points, sizeof(point_t) * num_points);
  lua_pushpointlist(lua, concat_data, len);
  return 1;
}

static int pointlist_len(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_POINT_LIST);
  lua_pushinteger(lua, var->size);
  return 1;
}

#if 0
static int pointlist_index(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_POINT_LIST);
  int i = (int)lua_tonumber(lua, -2);
  if (i < 0)
    i += var->size;
  point_t* data = var->datum;
  point_t* x = point_new(data[i].x, data[i].y, data[i].z);
  lua_pushpoint(lua, x);
  return 1;
}
#endif

void lua_pushpointlist(struct lua_State* lua, point_t* points, int size)
{
  // Bundle it up and store it in the given variable.
  store_pointlist(lua, points, size);
  lua_meta_key_val_t metatable[] = 
    {{"__tostring", pointlist_tostring},
     {"__concat", pointlist_concat},
     {"__len", pointlist_len},
//     {"__index", pointlist_index}, // FIXME: Currently broken!
     {NULL, NULL}};
  set_metatable(lua, "pointlist_metatable", metatable);
}

bool lua_isvector(struct lua_State* lua, int index)
{
  // A sequence with 3 numbers in it will work.
  if (lua_issequence(lua, index) && (lua_rawlen(lua, index) == 3))
    return true;
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_VECTOR);
}

vector_t* lua_tovector(struct lua_State* lua, int index)
{
  if (!lua_isvector(lua, index))
    return NULL;
  if (lua_istable(lua, index))
  {
    int len;
    real_t* seq = lua_tosequence(lua, index, &len);
    vector_t* v = vector_new(seq[0], seq[1], seq[2]);
    free(seq);
    return v;
  }
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_VECTOR)
    return (vector_t*)storage->datum;
  else if (storage->type == INTERPRETER_SEQUENCE)
  {
    real_t* seq = storage->datum;
    return vector_new(seq[0], seq[1], seq[2]);
  }
  else
    return NULL;
  if (!lua_isvector(lua, index))
    return NULL;
}

static int vector_tostring(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_VECTOR);
  vector_t* data = var->datum;
  char str[80];
  sprintf(str, "{%g, %g, %g}", data->x, data->y, data->z);
  lua_pushstring(lua, str);
  return 1;
}

static int vector_index(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_VECTOR);
  if (!lua_isnumber(lua, -1))
    return luaL_error(lua, "Non-numeric index given for vector.");
  int index = (int)lua_tonumber(lua, -1);
  if ((index < 1) || (index > 3))
    return luaL_error(lua, "Invalid index for vector: %d", index);
  vector_t* data = var->datum;
  real_t comp = (index == 1) ? data->x : (index == 2) ? data->y : data->z;
  lua_pushnumber(lua, comp);
  return 1;
}

static int vector_mul(lua_State* lua)
{
  interpreter_storage_t* var;
  int factor;
  if (lua_isnumber(lua, -1))
  {
    factor = (int)lua_tonumber(lua, -1);
    var = (void*)lua_topointer(lua, -2);
  }
  else
  {
    factor = (int)lua_tonumber(lua, -2);
    var = (void*)lua_topointer(lua, -1);
  }
  ASSERT(var->type == INTERPRETER_VECTOR);

  vector_t* data = var->datum;
  vector_t* veclist = malloc(sizeof(vector_t) * factor);
  for (int i = 0; i < factor; ++i)
    veclist[i] = *data;
  lua_pushvectorlist(lua, veclist, factor);
  return 1;
}

void lua_pushvector(struct lua_State* lua, vector_t* vec)
{
  // Bundle it up and store it in the given variable.
  store_vector(lua, vec);
  lua_meta_key_val_t metatable[] = 
    {{"__tostring", vector_tostring},
     {"__index", vector_index},
     {"__mul", vector_mul},
     {NULL, NULL}};
  set_metatable(lua, "vector_metatable", metatable);
}

bool lua_isvectorlist(struct lua_State* lua, int index)
{
  index = lua_absindex(lua, index);
  if (lua_istable(lua, index)) // A table with vectors in it will work.
  {
    size_t len = lua_rawlen(lua, index);
    if (len == 0)
      return false;
    for (size_t i = 1; i <= len; ++i)
    {
      lua_pushinteger(lua, (lua_Integer)i);
      lua_gettable(lua, index);
      bool is_vector = lua_isvector(lua, -1);
      lua_pop(lua, 1);
      if (!is_vector) 
        return false;
    }
    return true;
  }
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_POINT_LIST);
}

vector_t* lua_tovectorlist(struct lua_State* lua, int index, int* size)
{
  if (!lua_isvectorlist(lua, index))
    return NULL;
  index = lua_absindex(lua, index);
  if (lua_istable(lua, index))
  {
    *size = (int)lua_rawlen(lua, index);
    vector_t* vectors = malloc(sizeof(vector_t) * (*size));
    for (int i = 0; i < *size; ++i)
    {
      lua_pushinteger(lua, (lua_Integer)(i+1));
      lua_gettable(lua, index);
      vector_t* v = lua_tovector(lua, -1);
      lua_pop(lua, 1);
      vector_copy(&vectors[i], v);
    }
    return vectors;
  }
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_POINT_LIST)
  {
    *size = storage->size;
    return (vector_t*)storage->datum;
  }
  else
    return NULL;
}

#if 0
static int vectorlist_index(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_VECTOR_LIST);
  int i = (int)lua_tonumber(lua, -2);
  if (i < 0)
    i += var->size;
  vector_t* data = var->datum;
  vector_t* v = vector_new(data[i].x, data[i].y, data[i].z);
  lua_pushvector(lua, v);
  return 1;
}
#endif

void lua_pushvectorlist(struct lua_State* lua, vector_t* vectors, int size)
{
  // Bundle it up and store it in the given variable.
  store_vectorlist(lua, vectors, size);
  // Since point lists and vector lists are stored in an identical fashion, 
  // we can simply piggyback on the point list's metatable.
  lua_meta_key_val_t metatable[] = 
    {{"__tostring", pointlist_tostring},
     {"__concat", pointlist_concat},
     {"__len", pointlist_len},
//     {"__index", vectorlist_index}, // Slightly different! FIXME: Currently broken
     {NULL, NULL}};
  set_metatable(lua, "vectorlist_metatable", metatable);
}

bool lua_isboundingbox(struct lua_State* lua, int index)
{
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_BOUNDING_BOX);
}

bbox_t* lua_toboundingbox(struct lua_State* lua, int index)
{
  if (!lua_isuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_BOUNDING_BOX)
    return (bbox_t*)storage->datum;
  else
    return NULL;
}

static int boundingbox_tostring(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_BOUNDING_BOX);
  bbox_t* data = var->datum;
  char str[80];
  sprintf(str, "bounding_box{x1 = %g, x2 = %g, y1 = %g, y2 = %g, z1 = %g, z2 = %g}", 
    data->x1, data->x2, data->y1, data->y2, data->z1, data->z2);
  lua_pushstring(lua, str);
  return 1;
}

void lua_pushboundingbox(struct lua_State* lua, bbox_t* bbox)
{
  // Bundle it up and store it in the given variable.
  store_boundingbox(lua, bbox);
  lua_meta_key_val_t metatable[] = 
    {{"__tostring", boundingbox_tostring},
     {NULL, NULL}};
  set_metatable(lua, "boundingbox_metatable", metatable);
}

bool lua_isscalarfunction(struct lua_State* lua, int index)
{
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_SCALAR_FUNCTION);
}

st_func_t* lua_toscalarfunction(struct lua_State* lua, int index)
{
  if (!lua_isuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_SCALAR_FUNCTION)
    return (st_func_t*)storage->datum;
  else
    return NULL;
}

static int scalarfunction_tostring(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_SCALAR_FUNCTION);
  st_func_t* data = var->datum;
  char str[1024];
  sprintf(str, "scalar function (%s)", st_func_name(data));
  lua_pushstring(lua, str);
  return 1;
}

static int scalarfunction_call(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_SCALAR_FUNCTION);
  st_func_t* f = var->datum;
  ASSERT(st_func_num_comp(f) == 1);
  int num_args = lua_gettop(lua);
  if ((num_args != 2) && (num_args != 3))
    return luaL_error(lua, "Invalid argument(s). A scalar function takes x, t as arguments.");

  if (!lua_ispoint(lua, 2) && !lua_ispointlist(lua, 2))
    return luaL_error(lua, "Argument 1 must be a point or list of points.");

  if ((num_args == 3) && !lua_isnumber(lua, 3))
    return luaL_error(lua, "Argument 2, if given, must be a time.");

  if (lua_ispoint(lua, 2))
  {
    point_t* x = lua_topoint(lua, 2);
    real_t t = (real_t)lua_tonumber(lua, 3);
    real_t v;
    st_func_eval(f, x, t, &v);
    lua_pushnumber(lua, v);
  }
  else
  {
    int num_points;
    point_t* x = lua_topointlist(lua, 2, &num_points);
    real_t t = (real_t)lua_tonumber(lua, 3);
    real_t* v = malloc(sizeof(real_t) * num_points);
    for (int i = 0; i < num_points; ++i)
      st_func_eval(f, &x[i], t, &v[i]);
    lua_pushsequence(lua, v, num_points);
  }

  return 1;
}

void lua_pushscalarfunction(struct lua_State* lua, st_func_t* func)
{
  // Only single-component functions are allowed.
  ASSERT(st_func_num_comp(func) == 1); 
  // Bundle it up and store it in the given variable.
  store_scalar_function(lua, func);
  lua_meta_key_val_t metatable[] = 
    {{"__tostring", scalarfunction_tostring},
     {"__call", scalarfunction_call},
     {NULL, NULL}};
  set_metatable(lua, "scalarfunction_metatable", metatable);
}

bool lua_isvectorfunction(struct lua_State* lua, int index)
{
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_VECTOR_FUNCTION);
}

st_func_t* lua_tovectorfunction(struct lua_State* lua, int index)
{
  if (!lua_isuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_VECTOR_FUNCTION)
    return (st_func_t*)storage->datum;
  else
    return NULL;
}

static int vectorfunction_tostring(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_VECTOR_FUNCTION);
  st_func_t* data = var->datum;
  char str[1024];
  sprintf(str, "vector function (%s)", st_func_name(data));
  lua_pushstring(lua, str);
  return 1;
}

static int vectorfunction_call(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_VECTOR_FUNCTION);
  st_func_t* f = var->datum;
  ASSERT(st_func_num_comp(f) == 3);
  int num_args = lua_gettop(lua);
  if ((num_args != 2) && (num_args != 3))
    return luaL_error(lua, "Invalid argument(s). A vector function takes x, t as arguments.");

  if (!lua_ispoint(lua, 2) && !lua_ispointlist(lua, 2))
    return luaL_error(lua, "Argument 1 must be a point or a list of points.");

  if ((num_args == 3) && !lua_isnumber(lua, 3))
    return luaL_error(lua, "Argument 2, if given must be a time.");

  if (lua_ispoint(lua, 2))
  {
    point_t* x = lua_topoint(lua, 2);
    real_t t = (real_t)lua_tonumber(lua, 3);
    real_t v[3];
    st_func_eval(f, x, t, v);
    vector_t V = {.x = v[0], .y = v[1], .z = v[2]};
    lua_pushvector(lua, &V);
  }
  else
  {
    int num_points;
    point_t* x = lua_topointlist(lua, 2, &num_points);
    real_t t = (real_t)lua_tonumber(lua, 3);
    vector_t* V = malloc(sizeof(vector_t) * num_points);
    for (int i = 0; i < num_points; ++i)
    {
      real_t v[3];
      st_func_eval(f, x, t, v);
      V[i].x = v[0], V[i].y = v[1], V[i].z = v[2];
    }
    lua_pushvectorlist(lua, V, num_points);
  }

  return 1;
}

void lua_pushvectorfunction(struct lua_State* lua, st_func_t* func)
{
  // Only 3-component functions are allowed.
  ASSERT(st_func_num_comp(func) == 3); 
  // Bundle it up and store it in the given variable.
  store_vector_function(lua, func);
  lua_meta_key_val_t metatable[] = 
    {{"__tostring", vectorfunction_tostring},
     {"__call", vectorfunction_call},
     {NULL, NULL}};
  set_metatable(lua, "vectorfunction_metatable", metatable);
}

bool lua_issymtensorfunction(struct lua_State* lua, int index)
{
  if (lua_isnumber(lua, index)) 
    return true; // Diagonal symmetric tensor.
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return ((storage->type == INTERPRETER_SYM_TENSOR_FUNCTION) || 
          (storage->type == INTERPRETER_NUMBER));
}

st_func_t* lua_tosymtensorfunction(struct lua_State* lua, int index)
{
  if (lua_isnumber(lua, index))
  {
    real_t F = (real_t)lua_tonumber(lua, index);
    real_t v[6] = {F, 0.0, 0.0, F, 0.0, F};
    return constant_st_func_new(6, v);
  }
  if (!lua_isuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_SYM_TENSOR_FUNCTION)
    return (st_func_t*)storage->datum;
  else if (storage->type == INTERPRETER_NUMBER)
  {
    real_t F = *(real_t*)storage->datum;
    real_t v[6] = {F, 0.0, 0.0, F, 0.0, F};
    return constant_st_func_new(6, v);
  }
  else
    return NULL;
}

static int symtensorfunction_tostring(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_SYM_TENSOR_FUNCTION);
  st_func_t* data = var->datum;
  char str[1024];
  sprintf(str, "sym tensor function (%s)", st_func_name(data));
  lua_pushstring(lua, str);
  return 1;
}

static int symtensorfunction_call(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_SYM_TENSOR_FUNCTION);
  st_func_t* f = var->datum;
  ASSERT(st_func_num_comp(f) == 6);
  int num_args = lua_gettop(lua);
  if ((num_args != 2) && (num_args != 3))
    return luaL_error(lua, "Invalid argument(s). A sym tensor function takes x, t as arguments.");

  if (!lua_ispoint(lua, 2) && !lua_ispointlist(lua, 2))
    return luaL_error(lua, "Argument 1 must be a point or a list of points.");

  if ((num_args == 3) && !lua_isnumber(lua, 3))
    return luaL_error(lua, "Argument 2, if given must be a time.");

  if (lua_ispoint(lua, 2))
  {
    point_t* x = lua_topoint(lua, 2);
    real_t t = (real_t)lua_tonumber(lua, 3);
    real_t v[6];
    st_func_eval(f, x, t, v);
    lua_pushsequence(lua, v, 6);
  }
  else
  {
    int num_points;
    point_t* x = lua_topointlist(lua, 2, &num_points);
    real_t t = (real_t)lua_tonumber(lua, 3);
    real_t* V = malloc(sizeof(real_t) * 6 * num_points);
    for (int i = 0; i < num_points; ++i)
    {
      real_t v[6];
      st_func_eval(f, x, t, v);
      for (int j = 0; j < 6; ++j)
        V[6*i+j] = v[j];
    }
    lua_pushsequence(lua, V, 6*num_points);
  }

  return 1;
}

void lua_pushsymtensorfunction(struct lua_State* lua, st_func_t* func)
{
  // Only 6-component functions are allowed.
  ASSERT(st_func_num_comp(func) == 6); 
  // Bundle it up and store it in the given variable.
  store_sym_tensor_function(lua, func);
  lua_meta_key_val_t metatable[] = 
    {{"__tostring", symtensorfunction_tostring},
     {"__call", symtensorfunction_call},
     {NULL, NULL}};
  set_metatable(lua, "symtensorfunction_metatable", metatable);
}

bool lua_istensorfunction(struct lua_State* lua, int index)
{
  if (!lua_isnumber(lua, index))
    return true; // Diagonal symmetric tensor.
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return ((storage->type == INTERPRETER_TENSOR_FUNCTION) || 
          (storage->type == INTERPRETER_NUMBER));
}

st_func_t* lua_totensorfunction(struct lua_State* lua, int index)
{
  if (lua_isnumber(lua, index))
  {
    real_t F = (real_t)lua_tonumber(lua, index);
    real_t v[6] = {F, 0.0, 0.0, F, 0.0, F};
    return constant_st_func_new(6, v);
  }
  if (!lua_isuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_TENSOR_FUNCTION)
    return (st_func_t*)storage->datum;
  else if (storage->type == INTERPRETER_NUMBER)
  {
    real_t F = *(real_t*)storage->datum;
    real_t v[9] = {F, 0.0, 0.0, 0.0, F, 0.0, 0.0, 0.0, F};
    return constant_st_func_new(9, v);
  }
  else
    return NULL;
}

static int tensorfunction_tostring(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_TENSOR_FUNCTION);
  st_func_t* data = var->datum;
  char str[1024];
  sprintf(str, "tensor function (%s)", st_func_name(data));
  lua_pushstring(lua, str);
  return 1;
}

static int tensorfunction_call(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_TENSOR_FUNCTION);
  st_func_t* f = var->datum;
  ASSERT(st_func_num_comp(f) == 9);
  int num_args = lua_gettop(lua);
  if ((num_args != 2) && (num_args != 3))
    return luaL_error(lua, "Invalid argument(s). A tensor function takes x, t as arguments.");

  if (!lua_ispoint(lua, 2) && !lua_ispointlist(lua, 2))
    return luaL_error(lua, "Argument 1 must be a point or a list of points.");

  if ((num_args == 3) && !lua_isnumber(lua, 3))
    return luaL_error(lua, "Argument 2, if given must be a time.");

  if (lua_ispoint(lua, 2))
  {
    point_t* x = lua_topoint(lua, 2);
    real_t t = (real_t)lua_tonumber(lua, 3);
    real_t v[9];
    st_func_eval(f, x, t, v);
    lua_pushsequence(lua, v, 9);
  }
  else
  {
    int num_points;
    point_t* x = lua_topointlist(lua, 2, &num_points);
    real_t t = (real_t)lua_tonumber(lua, 3);
    real_t* V = malloc(sizeof(real_t) * 9 * num_points);
    for (int i = 0; i < num_points; ++i)
    {
      real_t v[9];
      st_func_eval(f, x, t, v);
      for (int j = 0; j < 9; ++j)
        V[9*i+j] = v[j];
    }
    lua_pushsequence(lua, V, 9*num_points);
  }

  return 1;
}

void lua_pushtensorfunction(struct lua_State* lua, st_func_t* func)
{
  // Only 6-component functions are allowed.
  ASSERT(st_func_num_comp(func) == 6); 
  // Bundle it up and store it in the given variable.
  store_tensor_function(lua, func);
  lua_meta_key_val_t metatable[] = 
    {{"__tostring", tensorfunction_tostring},
     {"__call", tensorfunction_call},
     {NULL, NULL}};
  set_metatable(lua, "tensorfunction_metatable", metatable);
}

bool lua_ismesh(struct lua_State* lua, int index)
{
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_MESH);
}

mesh_t* lua_tomesh(struct lua_State* lua, int index)
{
  if (!lua_isuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_MESH)
    return (mesh_t*)storage->datum;
  else
    return NULL;
}

static int mesh_tostring(lua_State* lua)
{
  interpreter_storage_t* var = (void*)lua_topointer(lua, -1);
  ASSERT(var->type == INTERPRETER_MESH);
  mesh_t* data = var->datum;
  char str[256];
  sprintf(str, "mesh (%d cells, %d faces, %d edges, %d nodes)", 
    data->num_cells, data->num_faces, data->num_edges, data->num_nodes);
  lua_pushstring(lua, str);
  return 1;
}

void lua_pushmesh(struct lua_State* lua, mesh_t* mesh)
{
  // Bundle it up and store it in the given variable.
  store_mesh(lua, mesh);
  lua_meta_key_val_t metatable[] = 
    {{"__tostring", mesh_tostring},
     {NULL, NULL}};
  set_metatable(lua, "mesh_metatable", metatable);
}

bool lua_isuserdefined(struct lua_State* lua, int index)
{
  if (!lua_isuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_USER_DEFINED);
}

void* lua_touserdefined(struct lua_State* lua, int index)
{
  if (!lua_isuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_USER_DEFINED)
    return (void*)storage->datum;
  else
    return NULL;
}

void lua_pushuserdefined(struct lua_State* lua, void* userdefined, void (*dtor)(void*))
{
  // Bundle it up and store it in the given variable.
  store_user_defined(lua, userdefined, dtor);
}

