from pybindgen import *

import sys
sys.path.append(".")
from PBGutils import *

#-------------------------------------------------------------------------------
# Class to handle wrapping this module.
#-------------------------------------------------------------------------------
class MeshEditor:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):
        
        # Includes
        mod.add_include('"PolytopeTypes.hh"')

        # Namespaces
        polytope = mod.add_cpp_namespace("polytope")

        # Expose types
        self.MeshEditor2d = addObject(polytope, "MeshEditor2d")
        self.MeshEditor3d = addObject(polytope, "MeshEditor3d")

        self.objs = [self.MeshEditor2d, self.MeshEditor3d]

        return
    
    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
                       
        for (obj, dim) in ((self.MeshEditor2d, 2),
                           (self.MeshEditor3d, 3)):
            self.generateMeshEditorBindings(obj, dim)

        return

    #---------------------------------------------------------------------------
    # Bindings (MeshEditor)
    #---------------------------------------------------------------------------
    def generateMeshEditorBindings(self, x, ndim):
        
        # Object names
        Tessellation = "polytope::Tessellation%id" % ndim
        
        # Constructors
        x.add_constructor([refparam(Tessellation, "mesh")])

        # Methods
        x.add_method("deleteCells", None,
                     [constrefparam("vector_of_unsigned", "cellsToDelete")])
        x.add_method("deleteFaces", None,
                     [constrefparam("vector_of_unsigned", "facesToDelete")])
        x.add_method("deleteNodes", None,
                     [constrefparam("vector_of_unsigned", "nodesToDelete")])
        x.add_method("cleanEdges", None,
                     [param("const double", "edgeTol")])
        
        # Attributes
        
        return
