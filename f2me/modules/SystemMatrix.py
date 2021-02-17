#My library                                                             
from modules.matrices_dir import matrix_Ax                                      
from modules.matrices_dir import matrix_Ay                                      
from modules.matrices_dir import matrix_PCoeff                                  
from modules.matrices_dir import matrix_Symm                                    
from modules.matrices_dir import matrix_Tn                                      
from modules.matrices_dir import matrix_Podd                                    
from modules.matrices_dir import matrix_Peven                                   
from modules.matrices_dir import matrix_BlockA                                  
from modules.matrices_dir import matrix_L
import numpy as np
import dolfin as df

class SystemMatrix:
    def __init__(self,my_input_cls, my_mesh_cls):
        #Input
        self.moment_order=my_input_cls['moment_order']
        self.Kn = my_input_cls['Kn']
        self.chi = my_input_cls['chi']
        self.epsilon_w = my_input_cls['epsilon_w']
        self.nv= my_mesh_cls.nv

    def read_system(self):
        A_x_list = (matrix_Ax.AxMatrix(self.moment_order))
        A_y_list = (matrix_Ay.AyMatrix(self.moment_order))
        P_Coeff_list = (matrix_PCoeff.PCoeffMatrix(self.moment_order,self.Kn))
        Symmetrizer_list = (matrix_Symm.SymmMatrix(self.moment_order))
        T_n_list = (matrix_Tn.TnMatrix(self.moment_order,self.nv))
        p_odd_list = (matrix_Podd.PoddMatrix(self.moment_order))
        p_even_list = (matrix_Peven.PevenMatrix(self.moment_order))
        Block_A_list = (matrix_BlockA.BlockAMatrix(self.moment_order))
        L_matrix_list = (matrix_L.LMatrix(self.moment_order,self.chi,self.epsilon_w))

        self.A_x = np.array(A_x_list)
        self.A_y = np.array(A_y_list)
        self.P_Coeff = np.array(P_Coeff_list)
        self.Symmetrizer = np.array(Symmetrizer_list)
        self.T_n = np.array(T_n_list)
        self.p_odd = np.array(p_odd_list)
        self.p_even = np.array(p_even_list)
        self.Block_A = np.array(Block_A_list)
        self.L_matrix = np.array(L_matrix_list)
        #Extra calculations
        self.Block_A_trans = np.transpose(Block_A_list)
        self.L_inverse = np.linalg.inv(L_matrix_list)

    def compute_system(self):
        self.SA_x = np.dot(self.Symmetrizer,self.A_x)
        self.SA_y = np.dot(self.Symmetrizer,self.A_y)
        self.SP_Coeff = np.dot(self.Symmetrizer,self.P_Coeff)
        #Calculate matrices of boundary integral term
        self.BC1 = np.dot(np.dot(np.dot(np.transpose(self.T_n),\
            np.transpose(self.p_odd)),self.L_inverse),np.dot(self.p_odd,self.T_n))
        self.BC2 = np.dot(np.dot(np.dot(np.dot(np.dot(np.transpose(self.T_n),\
            np.transpose(self.p_even)),self.Block_A_trans),self.L_matrix),\
            self.Block_A),np.dot(self.p_even,self.T_n))
        self.BC1_rhs = np.dot(np.transpose(self.T_n),\
                np.dot(np.transpose(self.p_even),self.Block_A_trans))
        self.BC2_rhs = np.dot(np.transpose(self.T_n),\
                np.dot(np.transpose(self.p_odd),self.L_inverse))

    def convert_to_ufl_form(self):
        #call the member functions
        self.read_system()
        self.compute_system()
        #
        self.SA_x_ufl = df.as_matrix(self.SA_x)
        self.SA_y_ufl = df.as_matrix(self.SA_y)
        self.SP_Coeff_ufl = df.as_matrix(self.SP_Coeff)
        self.BC1_ufl = df.as_matrix(self.BC1)
        self.BC2_ufl = df.as_matrix(self.BC2)
        self.BC1_rhs_ufl = df.as_matrix(self.BC1_rhs)
        self.BC2_rhs_ufl = df.as_matrix(self.BC2_rhs)


