import dolfin as df
class Solver():
    def __init__(self,my_input_cls, my_function_space_cls, my_var_prob_cls):
        self.my_input_cls = my_input_cls
        self.my_var_prob_cls = my_var_prob_cls
        self.my_function_space_cls = my_function_space_cls
        self.u = my_function_space_cls.u

    def inbuilt_newton_solver(self):
        V = self.my_function_space_cls.V
        F = self.my_var_prob_cls.F
        self.u = self.my_var_prob_cls.u
        du = df.TrialFunction(V)  #TrialFunctions --> wrong, without 's'
        Jac = df.derivative(F, self.u, du)

        #user given solver parameters
        abs_tol = self.my_input_cls['newton_abs_tol']
        rel_tol = self.my_input_cls['newton_rel_tol']
        max_itr = self.my_input_cls['newton_max_itr']
        step_size = self.my_input_cls['newton_relaxation_parameter']

        problem = df.NonlinearVariationalProblem(F,self.u,[],Jac)
        solver = df.NonlinearVariationalSolver(problem)
        solver.parameters ['newton_solver']['linear_solver'] = 'mumps'
        solver.parameters ['newton_solver']['absolute_tolerance'] = abs_tol
        solver.parameters ['newton_solver']['relative_tolerance'] = rel_tol
        solver.parameters ['newton_solver']['maximum_iterations'] = max_itr
        solver.parameters ['newton_solver']['relaxation_parameter'] =step_size
        #df.solve(F == 0, self.u, [], J=Jac,  solver_parameters={'newton_solver' : {'linear_solver' : 'mumps'}})
        solver.solve()
        #df.solve((a-L)==0, u,[])                                                   

    def custom_newton_solver(self):
        V = self.my_function_space_cls.V
        F = self.my_var_prob_cls.F
        #self.u = my_function_space_cls.u
        du = df.TrialFunction(V)  #TrialFunctions --> wrong, without 's'
        Jac = df.derivative(F, self.u, du)

        absolute_tol = 1E-5
        relative_tol = 9E-2
        u_inc = df.Function(V)
        nIter = 0
        relative_residual = 1
        CONVERGED =True
        MAX_ITER = 5000
        while relative_residual > relative_tol and nIter < MAX_ITER and CONVERGED:
            nIter += 1
            A, b = df.assemble_system(Jac, -(F), [])
            df.solve(A,u_inc.vector(),b)
            #df.solve(F == 0, u_inc.vector(), [], J=Jac,  solver_parameters={'newton_solver' : {'linear_solver' : 'mumps'}})
            relative_residual = u_inc.vector().norm('l2')
            #a = df.assemble(F)
            residual = b.norm('l2')
            lmbda = 0.8
            self.u.vector()[:] += lmbda*u_inc.vector()
            if residual > 1000:
                CONVERGED = False
                print("Did not converge after {} Iterations".format(nIter))
            elif residual < absolute_tol:
                CONVERGED =False
                print("converged after {} Iterations".format(nIter))         
            print ('{0:2d}  {1:3.2E}  {2:5e}'.format(nIter, relative_residual, residual))

    def inbuilt_linear_solver(self):
        a = self.my_var_prob_cls.a
        L = self.my_var_prob_cls.L
        V = self.my_function_space_cls.V
        self.u_Function = df.Function(V)
        df.solve(a==L, self.u_Function,[],solver_parameters={'linear_solver':'mumps'})

