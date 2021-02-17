import dolfin as df
import numpy
class FormulateVariationalProblem():
    def __init__(self,
            my_input_cls,
            my_mesh_cls,
            my_system_matrix_cls,
            my_function_space_cls
            ):
        #Input
        self.my_input_cls = my_input_cls
        self.my_system_matrix_cls = my_system_matrix_cls
        self.my_mesh_cls = my_mesh_cls
        self.my_function_space_cls = my_function_space_cls
        #Class Input
        self.moment_order = self.my_input_cls['moment_order']
        self.number_of_moments = self.my_input_cls['number_of_moments']
        self.problem_type = self.my_input_cls['problem_type']
        self.Kn = self.my_input_cls['Kn']
        self.Ma = self.my_input_cls['Ma']
        self.stab_enable = self.my_input_cls['stabilization']['enable']
        self.stab_type = self.my_input_cls['stabilization']['stab_type']
        self.DELTA_T = self.my_input_cls['stabilization']['cip']['DELTA_T']
        self.DELTA_P = self.my_input_cls['stabilization']['cip']['DELTA_P']
        self.DELTA_U = self.my_input_cls['stabilization']['cip']['DELTA_U']
        self.ht = self.my_input_cls['stabilization']['cip']['ht']
        self.hp = self.my_input_cls['stabilization']['cip']['hp']
        self.hu = self.my_input_cls['stabilization']['cip']['hu']

        self.epsilon_w = self.my_input_cls['epsilon_w']
        self.chi = self.my_input_cls['chi']
        self.theta_w_inner = self.my_input_cls['bc'][3000]['theta_w']
        self.u_t_w_inner = self.my_input_cls['bc'][3000]['u_t_w']
        self.u_n_w_inner = self.my_input_cls['bc'][3000]['u_n_w']
        self.p_w_inner = self.my_input_cls['bc'][3000]['p_w']
        self.theta_w_outer = self.my_input_cls['bc'][3100]['theta_w']
        self.u_t_w_outer = self.my_input_cls['bc'][3100]['u_t_w']
        self.u_n_w_outer = self.my_input_cls['bc'][3100]['u_n_w']
        self.p_w_outer = self.my_input_cls['bc'][3100]['p_w']

        #Class SystemMatrix 
        self.SA_x = self.my_system_matrix_cls.SA_x_ufl
        self.SA_y = self.my_system_matrix_cls.SA_y_ufl
        self.SP_Coeff = self.my_system_matrix_cls.SP_Coeff_ufl
        self.BC1 = self.my_system_matrix_cls.BC1_ufl
        self.BC2 = self.my_system_matrix_cls.BC2_ufl
        self.BC1_rhs = self.my_system_matrix_cls.BC1_rhs_ufl
        self.BC2_rhs = self.my_system_matrix_cls.BC2_rhs_ufl
        self.Symmetrizer= self.my_system_matrix_cls.Symmetrizer
        #Class FunctionSpace
        self.u = self.my_function_space_cls.u
        self.v = self.my_function_space_cls.v
        self.V = self.my_function_space_cls.V
        #Class Mesh
        self.mesh = self.my_mesh_cls.mesh
        self.nv = self.my_mesh_cls.nv
        self.ds = self.my_mesh_cls.ds
        self.dS = self.my_mesh_cls.dS

        #Output
        self.a = 0.0                                                                 
        self.L = 0.0
        self.F = 0.0

    def create_a(self):
        local_a =0.0
        df.ds = self.ds #Import <----
        local_a += (+0.5 * ((df.inner(self.v, (self.SA_x * df.Dx(self.u, 0)))) +\
                (df.inner(self.v, (self.SA_y * df.Dx(self.u, 1))))) ) * df.dx 
        local_a += (-0.5 * ((df.inner(df.Dx(self.v, 0),(self.SA_x * self.u))) +\
                (df.inner(df.Dx(self.v, 1),(self.SA_y * self.u)))) ) * df.dx 
        local_a += (+0.5 * df.inner(self.v, ((self.BC1 + self.BC2) * self.u)) ) * df.ds 
        local_a += (df.inner(self.v, (self.SP_Coeff * self.u)) ) * df.dx 
        return local_a
    def apply_stabilization(self):
        df.dS = self.dS #Important <----
        h_msh = df.CellDiameter(self.mesh)                                              
        h_avg = ( h_msh("+") + h_msh("-") )/2.0                                 
        #Output
        local_stab = 0.0
        if self.stab_type == 'cip':
            if self.moment_order==3 or self.moment_order=='ns':                               
                local_stab += self.DELTA_T * h_avg**self.ht * df.jump(df.grad(self.u[0]),self.nv)\
                        * df.jump(df.grad(self.v[0]),self.nv) * df.dS #cip for temp
            if self.moment_order==6:                                                     
                local_stab += self.DELTA_P * h_avg**self.hp * df.jump(df.grad(self.u[0]),self.nv)\
                        * df.jump(df.grad(self.v[0]),self.nv) * df.dS #cip for pressure
                local_stab += self.DELTA_U * h_avg**self.hu * df.jump(df.grad(self.u[1]),self.nv)\
                        * df.jump(df.grad(self.v[1]),self.nv) * df.dS #cip for velocity_x
                local_stab += self.DELTA_U * h_avg**self.hu * df.jump(df.grad(self.u[2]),self.nv)\
                        * df.jump(df.grad(self.v[2]),self.nv) * df.dS #cip for velocity_y
            if self.moment_order == 13:                                                    
                local_stab += self.DELTA_T * h_avg**self.ht * df.jump(df.grad(self.u[3]),self.nv)\
                        * df.jump(df.grad(self.v[3]),self.nv) * df.dS #cip for temp
                local_stab += self.DELTA_P * h_avg**self.hp * df.jump(df.grad(self.u[0]),self.nv)\
                        * df.jump(df.grad(self.v[0]),self.nv) * df.dS #cip for pressure
                local_stab += self.DELTA_U * h_avg**self.hu * df.jump(df.grad(self.u[1]),self.nv)\
                        * df.jump(df.grad(self.v[1]),self.nv) * df.dS #cip for velocity_x
                local_stab += self.DELTA_U * h_avg**self.hu * df.jump(df.grad(self.u[2]),self.nv)\
                        * df.jump(df.grad(self.v[2]),self.nv) * df.dS #cip for velocity_y
            if self.moment_order == 'grad13':                                                    
                local_stab += self.DELTA_T * h_avg**self.ht * df.jump(df.grad(self.u[3]),self.nv)\
                        * df.jump(df.grad(self.v[6]),self.nv) * df.dS #cip for temp
                local_stab += self.DELTA_P * h_avg**self.hp * df.jump(df.grad(self.u[0]),self.nv)\
                        * df.jump(df.grad(self.v[0]),self.nv) * df.dS #cip for pressure
                local_stab += self.DELTA_U * h_avg**self.hu * df.jump(df.grad(self.u[1]),self.nv)\
                        * df.jump(df.grad(self.v[1]),self.nv) * df.dS #cip for velocity_x
                local_stab += self.DELTA_U * h_avg**self.hu * df.jump(df.grad(self.u[2]),self.nv)\
                        * df.jump(df.grad(self.v[2]),self.nv) * df.dS #cip for velocity_y
        elif self.stab_type == 'gls':
            local_stab = (0.0001* h_msh * df.inner(self.SA_x * df.Dx(self.u,0)\
                    + self.SA_y * df.Dx(self.u,1) + self.SP_Coeff * self.u, self.SA_x * df.Dx(self.v,0)\
                    + self.SA_y * df.Dx(self.v,1) + self.SP_Coeff * self.v ) * df.dx)
        return local_stab
    def inhomogeneity(self):
        phi_local = df.Expression("atan2(x[1],x[0])",degree=2)
        if self.moment_order == 3:
            G_rhs_inner_list = [-2.0 * self.theta_w_inner, 0.0]
            G_rhs_outer_list = [-2.0 * self.theta_w_outer, 0.0]
        elif self.moment_order =='grad13':
            #Input
            u_n_w_inner = df.Expression("{}".format(self.u_n_w_inner),degree=2,phi=phi_local)
            u_n_w_outer = df.Expression("{}".format(self.u_n_w_outer),degree=2,phi=phi_local)
            u_t_w_inner = df.Expression("{}".format(self.u_t_w_inner),degree=2,phi=phi_local)
            u_t_w_outer = df.Expression("{}".format(self.u_t_w_outer),degree=2,phi=phi_local)
            #*******
            G_rhs_inner_list = [-self.chi * u_n_w_inner,-self.chi * u_t_w_inner, 0.0]
            G_rhs_outer_list = [-self.chi * u_n_w_outer,-self.chi * u_t_w_outer, 0.0]
        elif self.moment_order == 'ns':
            G_rhs_inner_list = [-1.0 *self.epsilon_w * self.theta_w_inner + self.u_n_w_inner\
                    ,-1.0 * self.u_t_w_inner]
            G_rhs_outer_list = [-1.0 *self.epsilon_w * self.theta_w_outer + self.u_n_w_outer\
                    ,-1.0 * self.u_t_w_outer]
        elif self.moment_order == 6:
            #Input
            u_n_w_inner = df.Expression("{}".format(self.u_n_w_inner),degree=2,phi=phi_local)
            u_n_w_outer = df.Expression("{}".format(self.u_n_w_outer),degree=2,phi=phi_local)
            u_t_w_inner = df.Expression("{}".format(self.u_t_w_inner),degree=2,phi=phi_local)
            u_t_w_outer = df.Expression("{}".format(self.u_t_w_outer),degree=2,phi=phi_local)
            p_w_inner = df.Expression("{}".format(self.p_w_inner),degree=2,phi=phi_local)
            p_w_outer = df.Expression("{}".format(self.p_w_outer),degree=2,phi=phi_local)
            #*******
            G_rhs_inner_list = [-self.epsilon_w * self.chi * p_w_inner + u_n_w_inner,\
                    -self.chi * u_t_w_inner, 0.0, 0.0] # inner
            G_rhs_outer_list = [-self.epsilon_w * self.chi * p_w_outer + u_n_w_outer,\
                    -self.chi * u_t_w_outer, 0.0, 0.0] # outer
        elif self.moment_order == 13:
            #Input
            theta_w_inner = df.Expression("{}".format(self.theta_w_inner),degree=2,phi=phi_local)
            theta_w_outer = df.Expression("{}".format(self.theta_w_outer),degree=2,phi=phi_local)
            u_n_w_inner = df.Expression("{}".format(self.u_n_w_inner),degree=2,phi=phi_local)
            u_n_w_outer = df.Expression("{}".format(self.u_n_w_outer),degree=2,phi=phi_local)
            u_t_w_inner = df.Expression("{}".format(self.u_t_w_inner),degree=2,phi=phi_local)
            u_t_w_outer = df.Expression("{}".format(self.u_t_w_outer),degree=2,phi=phi_local)
            p_w_inner = df.Expression("{}".format(self.p_w_inner),degree=2,phi=phi_local)
            p_w_outer = df.Expression("{}".format(self.p_w_outer),degree=2,phi=phi_local)
            #*******
            G_rhs_inner_list = [
                    - self.epsilon_w * self.chi * p_w_inner + u_n_w_inner,
                    - self.chi * u_t_w_inner,
                    - 2 * self.chi * theta_w_inner,
                    + (2/5) * self.chi * theta_w_inner,
                    - (1/5) * self.chi * theta_w_inner,
                    + self.chi * u_t_w_inner
                    ] # inner
            G_rhs_outer_list = [
                    - self.epsilon_w * self.chi * p_w_outer + u_n_w_outer,
                    - self.chi * u_t_w_outer,
                    - 2 * self.chi * theta_w_outer,
                    + (2/5) * self.chi * theta_w_outer,
                    - (1/5) * self.chi * theta_w_outer,
                    + self.chi * u_t_w_outer
                    ] # outer
        G_rhs_inner = df.as_vector(G_rhs_inner_list) #UFL form conversion
        G_rhs_outer = df.as_vector(G_rhs_outer_list)
        return G_rhs_inner, G_rhs_outer

    def create_bc(self):
        #Input
        G_rhs_inner, G_rhs_outer = self.inhomogeneity()
        local_bc = 0.0
        df.ds = self.my_mesh_cls.ds #<---- Important
        local_bc += (-0.5 * df.inner(self.v,((self.BC1_rhs-self.BC2_rhs) \
                * G_rhs_inner))  ) * df.ds(3000) #RHS-2
        local_bc += (-0.5 * df.inner(self.v,((self.BC1_rhs-self.BC2_rhs)\
                * G_rhs_outer))  ) * df.ds(3100) #RHS-2
        return local_bc
    def add_nonlinearity(self):
        #Output
        local_NL =0.0 
        if self.moment_order=='ns':
            local_NL = numpy.array(
                    [
                        df.Constant(0.0),
                        df.Constant(0.0),
                        df.Constant(0.0),                                           
                        -((2.0/3.0) * self.u[1]**2 -(1.0/3.0) * self.u[2]**2), 
                        -(self.u[1] * self.u[2]),
                        -(-(1.0/3.0) * self.u[1]**2 + (2.0/3.0) * self.u[2]**2)
                        ]
                    )
            local_NL = -(1/self.Kn)* local_NL
            local_NL = numpy.dot(self.Symmetrizer,local_NL)
            local_NL = df.as_vector(local_NL)
            local_NL = df.inner(self.v,local_NL) * df.dx 
        elif self.moment_order == 'grad13':
            local_NL = numpy.array(
                    [
                        df.Constant(0.0),
                        df.Constant(0.0),
                        df.Constant(0.0),
                        (self.Ma/self.Kn) * ((2.0/3.0) * self.u[1]**2 - (1.0/3.0) * self.u[2]**2), 
                        (self.Ma/self.Kn) * (self.u[1] * self.u[2]),
                        (self.Ma/self.Kn) * (-(1.0/3.0) * self.u[1]**2 + (2.0/3.0) * self.u[2]**2),
                        df.Constant(0.0),
                        (self.Ma/self.Kn) * ( ( (10/9) * self.u[6] * self.u[1] )
                            +  ( (1/3) * (self.u[3]*self.u[1] + self.u[4] * self.u[2]) ) ),
                        (self.Ma/self.Kn) * ( ( (10/9) * self.u[6] * self.u[2] )
                            +  ( (1/3) * (self.u[4]*self.u[1] + self.u[5] * self.u[2]) ) )
                        ] 
                    )
            local_NL = numpy.dot(self.Symmetrizer,local_NL)
            local_NL = df.as_vector(local_NL)
            local_NL = df.inner(self.v,local_NL) * df.dx 
        return local_NL

    def source_term(self):
        if self.moment_order==3:
            F_rhs = [0.0]*self.number_of_moments
            F_rhs[0] = df.Expression("2.0 - 1.0 * pow(sqrt(pow(x[0],2)+pow(x[1],2)),2)",degree=2)
        elif self.moment_order== 'nono6':
            F_rhs = [0.0]*NoV
        elif self.moment_order== 'nono13':
            F_rhs = [0.0]*NoV
            R= Expression("sqrt(pow(x[0],2)+pow(x[1],2))",degree=2)
            F_rhs[0] = Expression("1.0 * (1.0 - (5.0*pow(R,2))/(18.0*pow(kn,2))) * cos(phi)",R=R,kn=Kn,phi=phi,degree=2)
            F_rhs[3] = Expression("1.0 * (1.0 - (5.0*pow(R,2))/(18.0*pow(kn,2))) * cos(phi)",R=R,kn=Kn,phi=phi,degree=2)
            F_rhs[0] = Expression("2.0 - 1.0 * pow(sqrt(pow(x[0],2)+pow(x[1],2)),2)",degree=2)
            F_rhs[3] = Expression("2.0 - 1.0 * pow(sqrt(pow(x[0],2)+pow(x[1],2)),2)",degree=2)
        if self.moment_order == 3:
            F_rhs = df.as_vector(F_rhs) # converting to ufl vector
            local_source_term = df.inner(self.v, F_rhs) * df.dx
        return local_source_term


    #Left Hand Side
    def create_lhs(self):
        self.a = self.create_a()
        if self.stab_enable == True:
            if self.stab_type=='cip':
                self.a += self.apply_stabilization()
            elif self.stab_type == 'gls':
                self.a += df.lhs(self.apply_stabilization())

    def create_rhs(self):
        self.L = self.create_bc()
        if self.moment_order == 3:
            self.L += self.source_term()
        if self.stab_type == 'gls':
            self.L += df.rhs(self.apply_stabilization())
        if self.problem_type == 'nonlinear':
            self.L += self.add_nonlinearity()
            self.F = self.a - self.L
