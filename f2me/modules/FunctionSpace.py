import dolfin as df
class FunctionSpace():
    """FunctionSpace."""

    def __init__(self, my_input_cls, my_mesh_cls):
        """__init__.

        :param my_input_cls:
        :param my_mesh_cls:
        """
        #Inputs
        self.mesh = my_mesh_cls.mesh
        self.number_of_moment = my_input_cls['number_of_moments']
        self.problem_type = my_input_cls['problem_type']
        #Returns
        self.u = 0 
        self.v = 0
        self.V = 0
    def set_function_space(self):
        """set_function_space."""
        self.V = df.VectorFunctionSpace(self.mesh,'P',1,dim=self.number_of_moment)
        #Set test function(s)
        v_list = df.TestFunctions(self.V) 
        #Convert to ufl form
        self.v = df.as_vector(v_list) 
        #set trial function(s)
        if self.problem_type == 'nonlinear':
            u_list = df.Function(self.V) 
        elif self.problem_type == 'linear':
            u_list = df.TrialFunctions(self.V)
        #Convert to ufl form
        self.u = df.as_vector(u_list) 
