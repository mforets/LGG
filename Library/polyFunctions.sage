# returns a polytope in the form Ax <= b
def polytopeFromHrep(A, b):
    ieqs_list = list();
    for i in range(A.nrows()):
        ieqs_list.append(list(-A.row(i)))  #change in sign, necessary since Polyhedron receives Ax+b>=0
        ieqs_list[i].insert(0,b[i])
    P = Polyhedron(ieqs = ieqs_list, base_ring = RDF)
    return P

# extract inequalities from polytope
def polytopeVtoH(P):
    # add all constraints defined in P
    m = len(P.Hrepresentation())
    n = len(vector( P.Hrepresentation()[0] ))-1
    b = vector(RR, m)
    A = matrix(RR, m, n)
    P_gen = P.Hrep_generator();
    i=0
    for pigen in P_gen:
        pi_vec = pigen.vector()
        A.set_row(i, -pi_vec[1:len(pi_vec)])
        b[i] = pi_vec[0]
        i+=1;
    return [A, b]

# extract inequalities from polytope
def exactProjection(P):
    # add all constraints defined in P
    m = len(P.Hrepresentation())
    n = len(vector( P.Hrepresentation()[0] ))-1
    b = vector(RR, m)
    A = matrix(RR, m, n)
    P_gen = P.Hrep_generator();
    i=0
    for pigen in P_gen:
        pi_vec = pigen.vector()
        A.set_row(i, -pi_vec[1:len(pi_vec)])
        b[i] = pi_vec[0]
        i+=1;
    return [A, b]


# generate infinity ball
def Binfty(center, radius):
    ndim = len(center)

    A = matrix(RR, 2*ndim, ndim); b = vector(RR,2*ndim)

    count = 0
    for i in range(ndim):
        diri = zero_vector(RR,ndim)
        diri[i] = 1

        # external bound
        A.set_row(count, diri)
        b[count] = center[i] + radius
        count += 1
        
        # internal bound
        A.set_row(count, -diri)
        b[count] = -(center[i] - radius)
        count += 1
    
    
    P = polytopeFromHrep(A, b)
    
    return P


#Compute support function of a convex polytope (given as (A,b), assuming: Ax <= b)
def supp_fun_polyhedron(d, P, showOutput=1): 
    
    s_LP = MixedIntegerLinearProgram(maximization=True, solver = "GLPK")
    x = s_LP.new_variable(integer=False, nonnegative=False)

    # objective function
    obj = sum(d[i]*x[i] for i in range(len(d)))
    s_LP.set_objective(obj)
    
    # extract the constraints from P
    m = len(P.Hrepresentation())
    n = len(vector( P.Hrepresentation()[0] ))-1
    b = vector(RR, m)
    A = matrix(RR, m, n)
    P_gen = P.Hrep_generator();
    i=0
    for pigen in P_gen:
        pi_vec = pigen.vector()
        A.set_row(i, -pi_vec[1:len(pi_vec)])
        b[i] = pi_vec[0]
        i+=1;
        
    s_LP.add_constraint(A * x <= b);    
        
    if (showOutput):
        print '**** Solve LP (using GLPK) ****'    
        s_LP.show()
    
    oval = s_LP.solve()
    xopt = s_LP.get_values(x);
    
    if (showOutput):
        print 'Objective Value:', oval
        for i, v in xopt.iteritems():
            print 'x_%s = %f' % (i, v)
        print '\n'
    return oval, xopt


#Compute support function at d of an ellipse input as x^tr*Q*x <= 1
def supp_fun_ellipse(d, Q, showOutput=1): 
    if (Q.is_singular()):
        print 'error: input matrix is not invertible'
        return
    return sqrt(d.inner_product((~Q)*d))

# receive a matrix B and a polyhedron U and output the polyhedron V = BU
def matTimesPol(B, U):
    vertices_U = U.vertices_list()
    vertices_V = [np.dot(B, vertices_U[i]) for i in range(len(vertices_U))]
    V = Polyhedron(vertices = vertices_V, base_ring=RDF)
    return V


# compute the radius of polytope P. it computed in the supremum norm.
def polyRadius(P):
    n = len(vector( P.Hrepresentation()[0] ))-1
    
    r = 0
    for i in range(n):
        # generate canonical direction
        d = zero_vector(RR,n)
        d[i] = 1        
        [oval, xopt] = supp_fun_polyhedron(d, P, showOutput=0)
        asf = abs(oval)
        if (asf >= r):
            r = asf;
            
        # change sign
        d[i] = -1
        [oval, xopt] = supp_fun_polyhedron(d, P, showOutput=0)
        asf = abs(oval)
        if (asf >= r):
            r = asf;
        
    
    return r


