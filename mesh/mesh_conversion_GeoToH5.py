
#!/usr/bin/env python3

# pylint: disable=invalid-name

"""
Converter from geo-format to a mesh in h5-format.
Usage:
    python3 mesh_conversion_GeoToH5.py ring.geo ring0.h5 "-setnumber p 0"
    python3 mesh_conversion_GeoToH5.py ring.geo ring1.h5 "-setnumber p 1"
"""
import os
import sys
import dolfin as df

# Constants
GMSH_PATH = "gmsh"
GEO_NAME = "ring"

geo_input_file = sys.argv[1]
h5_output_file = sys.argv[2]
gmsh_arguments = sys.argv[3] 
tmp_name = "tmp"
# Create msh-mesh with Gmsh
os.system(
        "{} {} -2 -o {}.msh {}".format(
            GMSH_PATH, gmsh_arguments, tmp_name, geo_input_file
            )
        )

# Convert msh-mesh to xml-mesh
os.system("dolfin-convert {0}.msh {0}.xml".format(tmp_name))

# Delete msh-mesh
os.remove("{}.msh".format(tmp_name))

# Read xml-mesh
mesh = df.Mesh("{}.xml".format(tmp_name))
subdomains = df.MeshFunction(
        "size_t", mesh, "{}_physical_region.xml".format(tmp_name))
boundaries = df.MeshFunction(
        "size_t", mesh, "{}_facet_region.xml".format(tmp_name))

# Delete xml-mesh
os.remove("{}.xml".format(tmp_name))
os.remove("{}_physical_region.xml".format(tmp_name))
os.remove("{}_facet_region.xml".format(tmp_name))
print("done successfully")
# Write h5-mesh
file = df.HDF5File(mesh.mpi_comm(), h5_output_file, "w")
file.write(mesh, "/mesh")
file.write(subdomains, "/subdomains")
file.write(boundaries, "/boundaries")
