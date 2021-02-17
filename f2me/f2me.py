import sys
import dolfin as df                                                    
import numpy as np                                                      
import os                                                               
import time as time_module 
import matplotlib.pyplot as plt
import csv                                                              
import yaml #version higher than 5.0

df.parameters['ghost_mode']= 'shared_facet' #for parallel program          
#My library      
from modules.SystemMatrix import SystemMatrix # <--System Matrix class
from modules.Mesh import H5Mesh  # <-- Mesh class
from modules.FunctionSpace  import FunctionSpace #<--Function Space Class
from modules.FormulateVariationalProblem import FormulateVariationalProblem
from modules.Solver import Solver

#Inputs                                  
if len(sys.argv) >1 :
    user_given_input_file= sys.argv[1]
else:
    print("Provide an input file")
    quit()

with open(user_given_input_file,'r') as input_file:
    my_input_cls = yaml.load(input_file, Loader=yaml.FullLoader)

#projection = False
#manual = False
#automatic = True

for current_mesh in range(len(my_input_cls["mesh_list"])):
    #my_current_mesh = my_input_cls['mesh_list'][0] 
    my_current_mesh = my_input_cls['mesh_list'][current_mesh] 
    #step-1 : Convert current_mesh_name into a char_list    
    #step-2: Select the 4th-to-last element in the char_list 
    #Step-3: Convert the 4th-to-last element to integer     
    #Now I can run any single mesh I want in whichever order
    mesh_number=0  #For saving files according to mesh number
    mesh_number = int(list(my_current_mesh)[-4])  
    print(mesh_number)


    #MESH
    my_mesh_cls = H5Mesh(my_current_mesh) #<-- Class instance is created
    #print(my_current_mesh)  
    print("Maximum value of the mesh: ", my_mesh_cls.mesh.hmax()) 
    #df.plot(my_mesh_cls.mesh)


    #FUNCTION SPACE
    problem_type = my_input_cls['problem_type']
    number_of_moments = my_input_cls['number_of_moments']
    my_function_space_cls = FunctionSpace(my_input_cls, my_mesh_cls) #<-- class instance
    my_function_space_cls.set_function_space()

    #SYSTEM MATRIX
    my_system_matrix_cls = SystemMatrix(my_input_cls, my_mesh_cls) # <-- class instance
    my_system_matrix_cls.convert_to_ufl_form()

    #VARIATIONAL PROBLEM
    my_var_prob_cls = FormulateVariationalProblem(
        my_input_cls,
        my_mesh_cls,
        my_system_matrix_cls,
        my_function_space_cls
        ) #<-- class instance
    my_var_prob_cls.create_lhs()
    my_var_prob_cls.create_rhs()

    #SOLVER
    my_solver_cls = Solver(my_input_cls, my_function_space_cls, my_var_prob_cls) #<-- class instance
    if problem_type == 'nonlinear':
        my_solver_cls.inbuilt_newton_solver()
        u = my_solver_cls.u
    else:
        my_solver_cls.inbuilt_linear_solver()
        u = my_solver_cls.u_Function

    #===============
    #Post-processing
    #===============
    sol = u.split()

    def write_func(field_name, variable_name, mesh_num):
        xdmffile_u = df.XDMFFile(df.MPI.comm_world,
                              'results_mathematica/{0}_{1}.xdmf'
                              .format(variable_name,mesh_num))
        xdmffile_u.write(field_name)
        xdmffile_u.close()
    for i in range(my_input_cls["number_of_moments"]):
        write_func(sol[i], i, mesh_number)
    print("Program terminated successfully")
    moment_order = my_input_cls["moment_order"]

    def ErrorCalculation(exact_sol,numerical_sol):
        #interpolate on the mesh
        es = df.interpolate(exact_sol,my_function_space_cls.V.sub(0).collapse())
        ns = df.interpolate(numerical_sol,my_function_space_cls.V.sub(0).collapse())
        #compute values at the vertex
        err_linf = np.max(np.abs(es.compute_vertex_values()- ns.compute_vertex_values()))
        #err_l2 = np.linalg.norm(es.compute_vertex_values()-ns.compute_vertex_values())
        err_l2 = df.errornorm(es,ns,"L2")#fenics inbuilt norm calculator
        max_l2= np.linalg.norm(es.compute_vertex_values())
        max_linf = np.max(np.abs(es.compute_vertex_values())) or 1
        normalised_err_l2= err_l2/max_linf
        normalised_err_linf = err_linf/max_linf
        return  normalised_err_l2, normalised_err_linf
    def PressureErrorCalculation(exact_sol,numerical_sol1,numerical_sol2):
        #interpolate on the mesh
        es = interpolate(exact_sol, my_function_space_cls.V.sub(0).collapse())
        ns1 = interpolate(numerical_sol1, my_function_space_cls.V.sub(0).collapse())
        ns2 = interpolate(numerical_sol2, my_function_space_cls.V.sub(0).collapse())
        #compute values at the vertex
        #err_l2 = np.linalg.norm(es.compute_vertex_values()-
        #        (ns1.compute_vertex_values()+ns2.compute_vertex_values()))
        err_l2 = df.errornorm(es,ns1+ns2,"L2")#fenics inbuilt norm calculator
        err_linf = np.max(np.abs(es.compute_vertex_values()-
            (ns1.compute_vertex_values()+ns2.compute_vertex_values())))
        max_l2= np.linalg.norm(es.compute_vertex_values())
        max_linf = np.max(np.abs(es.compute_vertex_values())) or 1
        normalised_err_l2= err_l2/max_linf
        normalised_err_linf = err_linf/max_linf
        return  normalised_err_l2, normalised_err_linf

    if my_input_cls['moment_order']== 3:
            with open("01_coeffs.cpp", "r") as file:
                exact_solution_cpp_code = file.read()
            load_value = df.compile_cpp_code(exact_solution_cpp_code)
            #Temperature
            t_e = df.CompiledExpression(load_value.Temperature(),degree=2)
            t_l2,t_linf = ErrorCalculation(t_e,sol[0])
            #Heat flux
            sx_e = df.CompiledExpression(load_value.Heatfluxx(),degree=2)
            sx_l2,sx_linf = ErrorCalculation(sx_e,sol[1])
            sy_e = df.CompiledExpression(load_value.Heatfluxy(),degree=2)
            sy_l2,sy_linf = ErrorCalculation(sy_e,sol[2])
            errors =[
                    t_l2, t_linf,
                    sx_l2, sx_linf,
                    sy_l2,sy_linf ,
                     ]
            print(errors) 
    if my_input_cls['moment_order']== 'nono-6':
            with open("01_coeffs.cpp", "r") as file:
                exact_solution_cpp_code = file.read()
            load_value = compile_cpp_code(exact_solution_cpp_code)
            #Pressure
            #p_e = CompiledExpression(load_value.Pressure(),degree=2)
            #p_l2, p_linf = PressureErrorCalculation(p_e,sol[0],sol[3])
            #velocity
            ux_e = CompiledExpression(load_value.Velocityx(),degree=2)
            ux_l2,ux_linf = ErrorCalculation(ux_e,sol[1])
            uy_e = CompiledExpression(load_value.Velocityy(),degree=2)
            uy_l2,uy_linf = ErrorCalculation(uy_e,sol[2])
            #Pressure
            t_e = CompiledExpression(load_value.Pressure(),degree=2)
            t_l2,t_linf = ErrorCalculation(t_e,sol[0])
            #Stress
            sxx_e = CompiledExpression(load_value.Stressxx(),degree=2)
            sxx_l2,sxx_linf = ErrorCalculation(sxx_e,sol[4])
            sxy_e = CompiledExpression(load_value.Stressxy(),degree=2)
            sxy_l2,sxy_linf = ErrorCalculation(sxy_e,sol[5])
            syy_e = CompiledExpression(load_value.Stressyy(),degree=2)
            syy_l2,syy_linf = ErrorCalculation(syy_e,sol[6])
            #Heat flux
            #sx_e = CompiledExpression(load_value.Heatfluxx(),degree=2)
            #sx_l2,sx_linf = ErrorCalculation(sx_e,sol[7])
            #sy_e = CompiledExpression(load_value.Heatfluxy(),degree=2)
            #sy_l2,sy_linf = ErrorCalculation(sy_e,sol[8])
            errors =[
                    t_l2, t_linf,
                    ux_l2, ux_linf,
                     uy_l2,uy_linf ,
                     sxx_l2,sxx_linf,
                     sxy_l2,sxy_linf,
                     syy_l2,syy_linf
                     ]
            print(errors) 


#%%%%%%%%%%%%%%%%% Heat system Debug norm calculatation 
###########################################################
'''
    #interpolating the function values on function space
    def interpolate_func(exact_sol,numerical_sol):
        field_exact = df.interpolate(exact_sol, my_function_space_cls.V.sub(0).collapse())
        field_numeric = df.interpolate(numerical_sol,my_function_space_cls.V.sub(0).collapse())
        return field_exact, field_numeric
    ##
    #Calculating L2 and Linf error
    def error_calc(field_exact,field_numerical):
        max_field_exact = np.max(np.abs(field_exact.compute_vertex_values())) or 1
        err_L2 = df.errornorm(field_exact,field_numerical,"L2")
        err_linf = np.max(np.abs(
            field_exact.compute_vertex_values()
            - field_numerical.compute_vertex_values()))
        return max_field_exact,err_L2,err_linf

    ##OLD CODE BLOCK TO READ EXACT SOLUTION FROM CPP FILES
    with open("01_coeffs.cpp", "r") as file:
        exact_solution_cpp_code = file.read()
    load_value = df.compile_cpp_code(exact_solution_cpp_code)
    pressure_exact = df.CompiledExpression(load_value.Temperature(),degree=2)
    #sx_exact = CompiledExpression(load_value.VelocityX(),degree=2)
    #sy_exact = CompiledExpression(load_value.VelocityY(),degree=2)
    #pressure
    field_e_p,field_p = interpolate_func(pressure_exact,sol[0])
    #write_func(field_e_p,'e_p')
    max_exact_p,p_L2,p_linf = error_calc(field_e_p,field_p) 
    print("Theta L2,Linf :",p_L2/max_exact_p,p_linf/max_exact_p)
'''

