import random

# Define a function for computing the matrix M
def Matrix_Construction(N, Fq):
    row_num = 0
    Index = N**2 + N + 1 # Index of Gamma0 in SL3
    Number_of_Rows = 4*N**2 + 5*N + 4 # Number of rows needed in matrix

    P = ProjectiveSpace(2, GF(N))
    Points = P.rational_points()
    Dict = P.rational_points_dictionary()
    
    M = matrix(Fq, Number_of_Rows, Index, sparse=True) 

    # Condition 3: f(x:y:0) = 0 
    for pt in Points:
        if pt[2] == 0:
            M[row_num, Dict[pt]] = 1
            row_num += 1   

    # Condition 4: SUM_(z in P2(Z/NZ)) f(x:y:z) = 0 

    # x = 0, y = 0: Note that (0:0:z) for all z are in (0:0:1) class, so add N-1 to that location
    M[row_num, Dict[P(0, 0, 1)]] = N - 1
    row_num += 1
    
    # x = 0, y = 1:
    L_4a = []
    for i in range(1, N):
        p_4a = [0, i, 1]
        L_4a.append(P(p_4a))
        
    for pt in L_4a:
        M[row_num, Dict[pt]] = 1

    row_num += 1

    # x = 1
    L_4b = []
    for i in range(0, N): 
        for j in range(1, N):
            p_4b = [1, i, j]
            L_4b.append(P(p_4b))
        
        for pt in L_4b:
            M[row_num, Dict[pt]] = 1

        L_4b = []
        row_num += 1

    # Condition 1a: f(x:y:z) = f(-x:y:z)
    for pt in Points:
        p_negxyz = [-pt[0], pt[1], pt[2]]
        p_1a = P(p_negxyz)
        if pt == p_1a:   # If (x:y:z) = (-x:y:z), ignore since no new info gained
            continue
        else:
            M[row_num, Dict[pt]] = 1
            M[row_num, Dict[p_1a]] = -1
            row_num += 1 

    # Condition 1b: f(x:y:z) = f(z:x:y)
    for pt in Points:
        p_zxy = [pt[2], pt[0], pt[1]]
        p_1b = P(p_zxy)
        if pt == p_1b:   # If (x:y:z) = (z:x:y), ignore since no new info gained
            continue
        else:
            M[row_num, Dict[pt]] = 1
            M[row_num, Dict[p_1b]] = -1
            row_num += 1

    # Condition 1c: f(x:y:z) = -f(y:x:z)
    for pt in Points:
        p_yxz = [pt[1], pt[0], pt[2]]
        p_1c = P(p_yxz)
        if pt == p_1c:
            M[row_num, Dict[pt]] = 2
            row_num += 1     
        else:
            M[row_num, Dict[pt]] = 1
            M[row_num, Dict[p_1c]] = 1
            row_num += 1     

    # Condition 2: f(x:y:z) + f(-y:x-y:z) + f(y-x:-x:z) = 0
    for pt in Points: 
        p_negy = [-pt[1], pt[0] - pt[1], pt[2]]
        p_ynegx = [pt[1] - pt[0], -pt[0], pt[2]]
        p_2a = P(p_negy)
        p_2b = P(p_ynegx)
        if pt == p_2a == p_2b:
            M[row_num, Dict[pt]] = 3
            row_num += 1
        else:
            M[row_num, Dict[pt]] = 1
            M[row_num, Dict[p_2a]] = 1
            M[row_num, Dict[p_2b]] = 1
            row_num += 1
    
    return M

# Define a function for computing the dimension of the cuspidal cohomology for a given prime level N and finite field Fq
def Cuspidal_Cohomology_Dimension(N, Fq):
    Index = N**2 + N + 1 # Index of Gamma0 in SL3
    M = Matrix_Construction(N, Fq)
    M_Transpose = M.T
    Rank = M_Transpose.rank(algorithm="linbox")
    Corank = Index - Rank

    return N, Corank

# Define the van Geemen, van der Kallen, Top, Verberkmoes algorithm for reducing modular symbols
def mod_abs_value(int, modulus):
    Zm = Integers(modulus)
    if lift(Zm(int)) <= modulus/2:
        return Zm(int)
    else:
        return ( lift(Zm(int)) - modulus )

def modform_reduction(mod_symbol):
    nonzero_det_symbols = []

    Q_T = mod_symbol.T
    Q_T_ref = Q_T.echelon_form()
    Q_cef = Q_T_ref.T
        
    m1 = Q_cef[(0, 0)]
    m2 = Q_cef[(1, 1)]
    m3 = Q_cef[(2, 2)]

    a = Q_cef[(1, 0)]
    b = Q_cef[(2, 0)]
    c = Q_cef[(2, 1)]

    # Case 1: m = m1
    if ( m1 != 1 and m1 != -1 ):
        m = m1
        x = vector(ZZ, [1, 0, 0])

    # Case 2: m = m2
    elif ( m1 == 1 or m1 == -1 ) and ( m2 != 1 and m2 != -1 ):
        m = m2
        if m1 == 1:
            x = vector(ZZ, [mod_abs_value(-a, m), 1, 0])
        elif m1 == -1:
            x = vector(ZZ, [mod_abs_value(a, m), 1, 0])

    # Case 3: m = m3
    elif ( m1 == 1 or m1 == -1 ) and ( m2 == 1 or m2 == -1 ) and ( m3 != 1 and m3 != -1 ):
        m = m3
        if m1 == 1 and m2 == 1:
            x = vector(ZZ, [mod_abs_value(a*c - b, m), mod_abs_value(-c, m), 1])
        elif m1 == 1 and m2 == -1:
            x = vector(ZZ, [mod_abs_value(-a*c - b, m), mod_abs_value(c, m), 1])
        elif m1 == -1 and m2 == 1:
            x = vector(ZZ, [mod_abs_value(-a*c + b, m), mod_abs_value(-c, m), 1])
        elif m1 == -1 and m2 == -1:
            x = vector(ZZ, [mod_abs_value(a*c + b, m), mod_abs_value(c, m), 1])
    
    # Build vector v
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]

    t1 = x1/m
    t2 = x2/m
    t3 = x3/m

    q1 = mod_symbol[0]
    q2 = mod_symbol[1]
    q3 = mod_symbol[2]

    v1 = vector(QQ, t1*q1)
    v2 = vector(QQ, t2*q2)
    v3 = vector(QQ, t3*q3)

    v = v1 + v2 + v3

    Q1 = matrix(ZZ, [v, q2, q3])
    Q2 = matrix(ZZ, [q1, v, q3])
    Q3 = matrix(ZZ, [q1, q2, v])

    for test_matrix in [Q1, Q2, Q3]:
        if test_matrix.det() == 0:
            continue
        else: 
            nonzero_det_symbols.append(test_matrix)

    return nonzero_det_symbols

# Define a function for computing the action of Hecke operators on the cuspidal cohomology for a given prime level N, finite field Fq, and list of Hecke primes l_list.
def Compute_Hecke_Operators(N, Fq, l):
    P = ProjectiveSpace(2, GF(N))
    Dict = P.rational_points_dictionary()

    M = Matrix_Construction(N, Fq)
    Kernel = M.right_kernel()
    Basis = Kernel.basis()
    dim_U = len(Basis)

    # Build B_i coset reps for E_l
    Reps_E = []

    for a in range(0, l):
        for b in range(0, l):
            v1 = [1, 0, 0]
            v2 = [0, 1, 0]
            v3 = [a*N, b, l]
            Bi = matrix([v1, v2, v3])
            Reps_E.append(Bi)

    for c in range(0, l):
        v1 = [1, 0, 0]
        v2 = [c*N, l, 0]
        v3 = [0, 0, 1]
        Bi = matrix([v1, v2, v3])
        Reps_E.append(Bi)

    v1 = [l, 0, 0]
    v2 = [0, 1, 0]
    v3 = [0, 0, 1]
    Bi = matrix([v1, v2, v3])
    Reps_E.append(Bi)
    
    # Build the initial (potential) spanning set
    Q_j_list = []

    x_list = random.sample(range(1, N), 2)
    y_list = random.sample(range(1, N), 2)
    z_list = random.sample(range(1, N), 1)

    for x_j in x_list:
        for y_j in y_list:
            for z_j in z_list:
                y_over_x = mod(Integer(y_j)/Integer(x_j), N)
                z_over_x = mod(Integer(z_j)/Integer(x_j), N)

                Q_j_xy = matrix(ZZ, [[1, 0, 0], [x_j, 1, 0], [y_j, 0, 1]])
                Q_j_yz = matrix(ZZ, [[1, 0, 0], [y_j, 1, 0], [z_j, 0, 1]])
                Q_j_zx = matrix(ZZ, [[1, 0, 0], [z_j, 1, 0], [x_j, 0, 1]])
                Q_j_1 = matrix(ZZ, [[1, 0, 0], [y_over_x, 1, 0], [z_over_x, 0, 1]])
                Q_j_list.append([Q_j_xy, Q_j_yz, Q_j_zx, Q_j_1])

    start_var = True
    num_spanning_symbs_increment = 0
    Basis_Matrix = matrix(Fq, dim_U, 0)
    Tf_Matrix = matrix(Fq, dim_U, 0)
    Span_Set_Loc = -1

    while start_var == True or Basis_Matrix.rank() != dim_U:
        start_var = False

        Tf_Matrix = Tf_Matrix.augment(matrix(Fq, dim_U, len(Q_j_list)))
        Basis_Matrix = Basis_Matrix.augment(matrix(Fq, dim_U, len(Q_j_list)))

        for Qxyz in Q_j_list:
            Bi_sum_list = vector(Fq, dim_U, sparse=True)
            Qj_loc = -1
            Span_Set_Loc += 1
            for Qj in Qxyz:
                Qj_loc = Qj_loc + 1
                Qxyz_sum_list = vector(Fq, dim_U, sparse=True)
                for Bi in Reps_E: 
                    Qij_sum_list = vector(Fq, dim_U, sparse=True)

                    R = Qj*Bi

                    uni_symbols_R = []
                    starting_list = modform_reduction(R)
                    working_symbols = list(starting_list)

                    while starting_list != []: 
                        for test_matrix in starting_list: 
                            if test_matrix.det() == 1:
                                uni_symbols_R.append(test_matrix)
                                working_symbols.remove(test_matrix)
                            elif test_matrix.det() == -1:
                                working_symbols.remove(test_matrix)
                                row_test_matrix = matrix(ZZ, 3, [-1, 0, 0, 0, 1, 0, 0, 0, 1])*test_matrix
                                uni_symbols_R.append(row_test_matrix)
                            else:
                                output = modform_reduction(test_matrix)
                                index_to_replace = working_symbols.index(test_matrix)
                                num_items_to_replace = 1
                                working_symbols[index_to_replace:index_to_replace + num_items_to_replace] = output
                                
                        starting_list = working_symbols

                    for Mij in uni_symbols_R:
                        proj_point = P(Mij[(0,0)], Mij[(1,0)], Mij[(2,0)])
                        function_output_list = vector(Fq, dim_U, sparse=True)

                        for i in range(0, dim_U):
                            function_output_list[i] = Basis[i][Dict[proj_point]]
                            Qij_sum_list[i] = Qij_sum_list[i] + function_output_list[i]
                    
                    for i in range(0, dim_U):
                        if Qj_loc < 3:
                            Qxyz_sum_list[i] = Qxyz_sum_list[i] + Qij_sum_list[i]
                        elif Qj_loc == 3:
                            Qxyz_sum_list[i] = Qxyz_sum_list[i] - Qij_sum_list[i]

                for j in range(0, dim_U):
                    Bi_sum_list[j] = Bi_sum_list[j] + Qxyz_sum_list[j]

            for k in range(0, dim_U):
                Tf_Matrix[k, Span_Set_Loc] = Bi_sum_list[k]
                Basis_Matrix[k, Span_Set_Loc] = Basis[k][Dict[P(list(Qxyz[0].column(0)))]] + Basis[k][Dict[P(list(Qxyz[1].column(0)))]] + Basis[k][Dict[P(list(Qxyz[2].column(0)))]] - Basis[k][Dict[P(list(Qxyz[3].column(0)))]]

        if Basis_Matrix.rank() != dim_U:
            num_spanning_symbs_increment += 1

            x_rand = choice(range(1, N))
            y_rand = choice(range(1, N))
            z_rand = choice(range(1, N))

            y_over_x = mod(Integer(y_rand)/Integer(x_rand), N)
            z_over_x = mod(Integer(z_rand)/Integer(x_rand), N)

            Q_j_xy = matrix(ZZ, [[1, 0, 0], [x_rand, 1, 0], [y_rand, 0, 1]])
            Q_j_yz = matrix(ZZ, [[1, 0, 0], [y_rand, 1, 0], [z_rand, 0, 1]])
            Q_j_zx = matrix(ZZ, [[1, 0, 0], [z_rand, 1, 0], [x_rand, 0, 1]])
            Q_j_1 = matrix(ZZ, [[1, 0, 0], [y_over_x, 1, 0], [z_over_x, 0, 1]])
            Q_j_list.append([Q_j_xy, Q_j_yz, Q_j_zx, Q_j_1])

        else:
            A = Basis_Matrix.solve_left(Tf_Matrix)
            Char_Poly = A.charpoly('T')
            A_Data = A.eigenspaces_right(format="galois")

            return l, Char_Poly, A_Data