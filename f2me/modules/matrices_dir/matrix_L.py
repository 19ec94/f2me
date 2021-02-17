#Store Podd matrices
def LMatrix(moment_order,chi,epsilon_w):
    if moment_order==3:
        L_Matrix = [
                [2.0*chi, 0.0    ],
                [0.0,     11.0*chi]
                ]
        return L_Matrix 
    elif moment_order=='grad13':
        L_Matrix = [
                [chi, 0.0, 0.0],
                [0.0, chi, 0.0],
                [0.0, 0.0, 2*chi]
                ]
        return L_Matrix 
    elif moment_order==6:
        L_Matrix = [
                [chi*epsilon_w, 0.0, 0.0, 0.0],
                [0.0, chi, 0.0, 0.0],
                [0.0, 0.0, (28.0/15.0)*chi, -(14.0/15.0)*chi],
                [0.0, 0.0, -(14.0/15.0)*chi, (22.0/15.0)*chi]
                ]
        return L_Matrix 
    elif moment_order==13:
        L_Matrix = [
                [epsilon_w*chi, 0, 0, 0, 0, 0],
                [0, chi, 0, 0, 0, -chi],
                [0, 0, 2*chi, -(2/5)*chi, (1/5)*chi, 0],
                [0, 0, -(2/5)*chi, (52/25)*chi, -(26/25)*chi, 0],
                [0, 0, +(1/5)*chi, -(26/25)*chi, (38/25)*chi, 0],
                [0, -chi,0, 0, 0, 13*chi]
                ]
        return L_Matrix
    elif moment_order=='ns':
        L_Matrix = [
                [epsilon_w*chi, 0.0    ],
                [0.0,     chi]
                ]
        return L_Matrix 
    elif moment_order=='g26w':
        L_Matrix = [
                [1, 0, 0.5, -0.4, 0.2, 0], 
                [0, 1, 0, 0, 0, 0.5], 
                [0.5, 0, 2.25, 0.2, -0.1, 0], 
                [-0.4, 0, 0.2, 2.24, -1.12, 0], 
                [0.2, 0, -0.1, -1.12, 1.56, 0], 
                [0, 0.5, 0, 0, 0, 3.25]
                ]
        return L_Matrix 
    else:
        print("Wrong L matrix")
