import dolfin as df 
class H5Mesh:
    """H5Mesh."""

    def __init__(self, h5_file):
        """__init__.

        :param h5_file:
        """
        self.mesh = df.Mesh()
        hdf =  df.HDF5File(self.mesh.mpi_comm(),h5_file,"r")
        hdf.read(self.mesh,"/mesh", False)
        dim = self.mesh.topology().dim()
        
        self.subdomains = df.MeshFunction("size_t", self.mesh, dim)
        hdf.read(self.subdomains, "/subdomains")
        
        self.boundaries = df.MeshFunction("size_t", self.mesh, dim-1)
        hdf.read(self.boundaries, "/boundaries")
        #Used in the program
        self.nv = df.FacetNormal(self.mesh) #normal vector
        self.ds = df.Measure("ds",domain=self.mesh, subdomain_data=self.boundaries)               
        self.dS = df.Measure("dS",domain=self.mesh, subdomain_data=self.boundaries) 
