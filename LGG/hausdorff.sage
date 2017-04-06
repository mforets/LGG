r"""
Extension to polyhedra library for the purpose of mathematical modeling, with a focus on computational geometry.

hausdorff.sage stands for Hausdorff distance approximation.
Contains functions to estimate Hausdorff distance between polytopes in different setups.

AUTHORS:

- Marcelo Forets (last modified 2016-10-20)

- Frederic Viry

EXAMPLES::

See 'Examples-polyFunctions-core.ipynb'
See '~/Projects/Compositional/dHp qblocks some tests.ipynb'

"""

#*****************************************************************************
#       Copyright (C) 2016 Marcelo Forets <mforets@nonlinearnotes.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# compatibility with older version (in SMC)
try:
    from polyFunctions_core_future import *
except ImportError:
    try:
        from polyFunctions_core import *
    except ImportError:
        raise ImportError('Could not import module polyFunctions_core.')

from norms import *

class NotImplementedException(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg


def asphericity_polytope(P, tol=1e-1, MaxIter=1e2, verbose=0, solver='GLPK'):
    r""" Algorithm for computing asphericity

    INPUT:

    "x0" - point in the interior of P

    "tol" - stop condition

    OUTPUT:

    "psi_opt" -  min_{x in P} psi(x) = R(x)/r(x), where R(x) and r(x) are the circumradius
                and the inradius at x respectively.

    "alphakList" - sequence of psi(x) which converges to psi_opt

    "xkList" - sequence of interior points which converge to the minimum

    """

    # check low-dimensionality
    if not P.is_full_dimensional():
        raise NotImplementedException('The polytope is not full-dimensional')

    # initial data
    p = P.ambient_dim()

    [DP, d] = PolyhedronToHSpaceRep(P) # in the form Dy <= d
    D = -DP    # transform to <Dj,y> + dj >= 0
    l = D.nrows()
    mP = opposite_polyhedron(P)

    # unit ball, infinity norm
    Bn = BoxInfty(center=zero_vector(RDF,p), radius=1)
    [C, c] = PolyhedronToHSpaceRep(Bn)
    C = -C    # transform to <Cj,y> + cj >= 0
    m = len(c)

    # auxiliary matrices A, a, B, b
    A = matrix(RDF, C.nrows(), C.ncols())
    a = []
    for i in range(m):
        A.set_row(i, -C.row(i)/c[i])
        a.append(supp_fun_polyhedron(mP, A.row(i),solver=solver))

    B = matrix(RDF, D.nrows(), D.ncols())
    b = []
    for j in range(l):
        [nDj, ymax] = supp_fun_polyhedron(Bn, D.row(j), return_xopt=True,solver=solver)
        B.set_row(j, D.row(j)/nDj)
        b.append(d[j]/nDj)

    # generate an initial point
    x0 = chebyshev_center(DP, d)

    # compute asphericity at x0
    R = max(A.row(i)*vector(x0)+a[i] for i in range(m))
    r = min(B.row(j)*vector(x0)+b[j] for j in range(l))
    alpha0 = R/r

    xkList = [vector(x0)]
    alphakList = [alpha0]

    xk0 = x0; psi_k0 = alpha0
    alphak0 = alpha0
    convergeFlag = 0
    iterCount = 0
    while (iterCount < MaxIter) and (not convergeFlag):
        [psi_k, xkDict, zkDict] = _asphericity_iteration(xk0, l, m, p, A, a, B, b, verbose=verbose,solver=solver)
        xk = vector(xkDict); xkList.append(xk)
        zk = vector(zkDict); alphakList.append(zk[0])
        xk0 = xk;
        if (abs(psi_k - psi_k0) < tol):
            convergeFlag = 1
        psi_k0 = psi_k
        iterCount += 1

    x_asph = xkList.pop()
    R_asph = max(A.row(i)*vector(x_asph)+a[i] for i in range(m))
    r_asph = min(B.row(j)*vector(x_asph)+b[j] for j in range(l))
    asph = R_asph/r_asph

    return [asph, x_asph]


def _asphericity_iteration(xk, l, m, p, A, a, B, b, verbose=0, solver='GLPK'):
    # what is the best way to share this data?

    #alphak = asphericity(xk)
    R_xk = max(A.row(i)*vector(xk)+a[i] for i in range(m))
    r_xk = min(B.row(j)*vector(xk)+b[j] for j in range(l))
    alphak = R_xk/r_xk

    a_LP = MixedIntegerLinearProgram(maximization=False, solver=solver)
    x = a_LP.new_variable(integer=False, nonnegative=False)
    z = a_LP.new_variable(integer=False, nonnegative=False)

    for i in range(m):
        for j in range(l):
            aux = sum( (A.row(i)[k]-alphak*B.row(j)[k])*x[k] for k in range(p))
            a_LP.add_constraint(aux + a[i] - alphak*b[j] <= z[0]);

    #solve LP
    a_LP.set_objective(z[0])

    if (verbose):
        print '**** Solve LP (using GLPK) ****'
        a_LP.show()

    opt_val = a_LP.solve()

    x_opt = a_LP.get_values(x);
    z_opt = a_LP.get_values(z)

    if (verbose):
        print 'Objective Value:', opt_val
        for i, v in x_opt.iteritems():
            print 'x_%s = %f' % (i, v)

        for i, v in z_opt.iteritems():
            print 'z = %f' % (v)
        print '\n'

    return [opt_val, x_opt, z_opt]


# obtain a sampling point from the interior of P
# check idea in cprnd.m
#def polyRandomPoint(P):
#    return x

# generate Gram-Schmidt orthogonalization process from columns of matrix A
def gram_schmidt_columns(A):
    import numpy as np

    Q, R = np.linalg.qr(A)
    return Q

# generate Gram-Schmidt orthogonalization process from rows of matrix A
def gram_schmidt_rows(A):
    import numpy as np

    Q, R = np.linalg.qr(A.T)
    return Q


def dHp(X, Xhat, p=2, verbose = 1):
    """ Compute approximation to Hausdorff distance using the rows of X, in the p-norm
    Use p=2 for Euclidean norm (default); p='inf' for sup-norm
    """
    [A, b] = polytopeToHrep(X)
    [C, d] = polytopeToHrep(Xhat)

    aux = 0;
    for i in range(A.nrows()):
        fi = A.row(i)
        sf = supp_fun_Ab(fi,C,d,verbose) #the optimum value is the first component, the second component is a dictionary containing the solution
        cfraction = (sf[0]-b[i])/vecpnorm(fi,p)
        if (cfraction > aux):
            aux = cfraction;

    pX = 1
    return aux*pX


#Compute dHp(X,Xhat) using the formula above. Assume that n = n1+n2. Assume that pX = 1
def dHp_2blocks(X, S, S1, S2, n1, n, p=2, verbose=1):

    [F, g] = polytopeToHrep(X)
    aux = 0;
    for i in range(F.nrows()):
        fi = F.row(i)

        BM1 = block_matrix(n1,n,[identity_matrix(n1),zero_matrix(n1,n2)])
        fi1 = S1.T*BM1*((~S).T*fi)
        sf1 = supp_fun_Ab(fi1,F,g,verbose)

        BM2 = block_matrix(n2,n,[zero_matrix(n2,n1), identity_matrix(n2)])
        fi2 = S2.T*BM2*((~S).T*fi)
        sf2 = supp_fun_Ab(fi2,F,g,verbose)

        cfraction = (sf1[0] + sf2[0] - g[i])/vecpnorm(fi,p)
        if (cfraction > aux):
            aux = cfraction;

        return aux


#Compute dHp(X,Xhat) using the formula above. nL = [n1,...,nq], with n = sum(nL). Assume that pX0 = 1
#By default we use p=2 norm; use p='inf' for sup-norm
def dHp_qblocks(X, S, nL, p=2, verbose=0):

    if 'list' in str(type(X)):
        F = X[0]; g = X[1]
    elif 'polyhedron' in str(type(X)):
        [F, g] = polytopeToHrep(X)


    if (verbose):
        print F # en el ejemplo, F tiene que ser de 12x6


    q = len(nL)
    n = sum(nL)

    aux = 0;
    for i in range(F.nrows()):
        fi = F.row(i)

        #print 'fi = ', fi

        aux_sfSum = 0
        for j in range(q):

            ILj = list()

            for k in [0..j-1]:
                ILj.append(zero_matrix(nL[j],nL[k]))

            ILj.append(identity_matrix(nL[j]))

            for k in [j+1..q-1]:
                ILj.append(zero_matrix(nL[j],nL[k]))

            I_nj_n = block_matrix(RR, ILj, nrows=1, ncols=q)
            #print 'Injn = ', I_nj_n

            #alternativo, mas corto, revisar:
            #I_nj_n = zero_matrix(RR, nL[j], n)
            #I_nj_n.set_block(0, nL[j], identity_matrix(nL[j]))


            Sj = S.submatrix(row = sum(nL[k] for k in [0..j-1]), col=0, nrows=nL[j])
            #print 'Sj = ', Sj

            #print 'before fij'
            #print len(fi), size(S)
            #type(I_nj_n)
            #type(S)
            #type(fi)

            fij = Sj.T*I_nj_n*((~S).T*fi)
            #print 'fij = ', fij

            #print 'after fij'
            sfj = supp_fun_Ab(fij, F, g, verbose)
            #print 'sfj = ', sfj

            aux_sfSum += sfj[0]


        cfraction = (aux_sfSum - g[i])/vecpnorm(fi,p)
        if (cfraction > aux):
            aux = cfraction;

    pX0 = 1

    return aux*pX0


# exploring automatically all partitions
# Warning: this search has exponential cost. By assumption we exclude q=1.
# By default we use p=2 norm; use p='inf' for sup-norm
def find_minimal_qBlocks(constraint_pairs_list, X, S, p=2, verbose = 0):

    if 'list' in str(type(X)):
        F = X[0]; g = X[1]
    elif 'polyhedron' in str(type(X)):
        [F, g] = polytopeToHrep(X)

    r = len(constraint_pairs_list)
    n = F.ncols()

    # generate sticks ...|..|...|..
    allowed_list = [1+i for i in range(n-1)]
    for i in range(r):
        allowed_list.remove(constraint_pairs_list[i][0])

    print 'allowed list (position of sticks) = ', allowed_list

    dHp_all = list()
    nL_all = list()

    # generate all combinations of q blocks
    for q in [2..n-r]: #q-1 is the number of sticks

        if (verbose):
            print '\n', 'q ', q

        C = Combinations(allowed_list, q-1)
        CList = C.list()

        for k in range(len(CList)): #k is an index for each possible combination
            nL_aux = list()
            nL_aux.append(CList[k][0])
            for j in [1..q-2]:
                nL_aux.append(CList[k][j]-CList[k][j-1])
            nL_aux.append(n - CList[k][q-2])

            dHp_aux = dHp_qblocks(X, S, nL_aux, p, verbose=0)
            nL_all.append(nL_aux)
            dHp_all.append(dHp_aux)

            if (verbose):
                print 'k ', k, '  Comb = ', CList[k],  '  nL = ', nL_aux, '  dHp = ', dHp_aux.N(digits=6)


    return dHp_all, nL_all



#Compute dHp(X,Xhat) using the formula above. nL = [n1,...,nq], with n = sum(nL). Assume that pX0 = 1
#By default we use p=2 norm; use p='inf' for sup-norm
def dHp_qblocks_optimized(F, g, S, nL, p=2, verbose=0):

    q = len(nL)
    n = sum(nL)
    #print 'q=', q
    #print 'n=', n

    nL.insert(0,0);
    indices_partitions = matrix(ZZ, np.tril(np.ones((q+1,q+1)), 0))*vector(ZZ, nL);

    # check that the size of S equals n
    # how to use Cython in these loops?
    # how to incorporate Gurobi solver in the supp_fun_Ab(..)?

    aux = 0;

    # Prepare LP
    s_LP = MixedIntegerLinearProgram(maximization=True, solver = "GLPK")
    x = s_LP.new_variable(integer=False, nonnegative=False)

    # numpy arrays : use a different solver !
    #s_LP.add_constraint(matrix(F) * x <= vector(g)); # without the cast to matrix/vector, I get an unsupported operand error

    # sage dense matrices
    s_LP.add_constraint(F * x <= g);


    for i in range(F.nrows()):
        fi = F.row(i)

        aux_sfSum = 0

        #print 'i=', i

        for j in range(q):

            #ILj = list()
            #for k in [0..j-1]:
            #    ILj.append(zero_matrix(nL[j],nL[k]))
            #ILj.append(identity_matrix(nL[j]))
            #for k in [j+1..q-1]:
            #    ILj.append(zero_matrix(nL[j],nL[k]))
            #I_nj_n = block_matrix(RR, ILj, nrows=1, ncols=q)
            #print I_nj_n
            #Sj = S.submatrix(row = sum(nL[k] for k in [0..j-1]), col=0, nrows=nL[j])
            #fij = Sj.T*I_nj_n*((~S).T*fi)

            # if S is a sage dense matrix
            M = S[indices_partitions[j]:indices_partitions[j+1],:].T*S[indices_partitions[j]:indices_partitions[j+1],:];
            fij = M*fi;

            # Set objective function and solve LP
            # sfj = supp_fun_Ab(fij, F, g, verbose)
            obj = sum(fij[i]*x[i] for i in range(n))
            s_LP.set_objective(obj)

            sfj = s_LP.solve()

            aux_sfSum += sfj

            #print '    j= ', j, '    sfj = ', sfj[0]
            #print '     fij = ', fij

        # define next term
        #cfraction = (aux_sfSum - g[i])/vecpnorm(fi,p)
        # Recall that the coeff matrix satisfies vecpnorm(fi, p)=1 for most models.
        cfraction = (aux_sfSum - g[i])

        if (cfraction > aux):
            aux = cfraction;

    #pX0 = 1 # asphericity is 'assumed' = 1
    # aux = aux*px0
    return aux


#Compute dHp(X,Xhat) using the formula above. nL = [n1,...,nq], with n = sum(nL). Assume that pX0 = 1
#By default we use p=2 norm; use p='inf' for sup-norm
def dHp_qblocks_optimized_numpy(F, g, S, nL, p=2, verbose=0):

    q = len(nL)
    n = sum(nL)
    #print 'q=', q
    #print 'n=', n

    nL.insert(0,0);
    indices_partitions = matrix(ZZ, np.tril(np.ones((q+1,q+1)), 0))*vector(ZZ, nL);

    # check that the size of S equals n
    # how to use Cython in these loops?
    # how to incorporate Gurobi solver in the supp_fun_Ab(..)?

    aux = 0;

    # Prepare LP
    s_LP = MixedIntegerLinearProgram(maximization=True, solver = "GLPK")
    x = s_LP.new_variable(integer=False, nonnegative=False)

    # numpy arrays : use a different solver !
    #s_LP.add_constraint(matrix(F) * x <= vector(g)); # without the cast to matrix/vector, I get an unsupported operand error

    # sage dense matrices
    s_LP.add_constraint(F * x <= g);


    for i in range(F.nrows()):
        fi = F.row(i)

        aux_sfSum = 0

        #print 'i=', i

        for j in range(q):

            #ILj = list()
            #for k in [0..j-1]:
            #    ILj.append(zero_matrix(nL[j],nL[k]))
            #ILj.append(identity_matrix(nL[j]))
            #for k in [j+1..q-1]:
            #    ILj.append(zero_matrix(nL[j],nL[k]))
            #I_nj_n = block_matrix(RR, ILj, nrows=1, ncols=q)
            #print I_nj_n
            #Sj = S.submatrix(row = sum(nL[k] for k in [0..j-1]), col=0, nrows=nL[j])
            #fij = Sj.T*I_nj_n*((~S).T*fi)

            # if S is a sage dense matrix
            M = S[indices_partitions[j]:indices_partitions[j+1],:].T*S[indices_partitions[j]:indices_partitions[j+1],:];
            fij = M*fi;

            # Set objective function and solve LP
            # sfj = supp_fun_Ab(fij, F, g, verbose)
            obj = sum(fij[i]*x[i] for i in range(n))
            s_LP.set_objective(obj)

            sfj = s_LP.solve()

            aux_sfSum += sfj

            #print '    j= ', j, '    sfj = ', sfj[0]
            #print '     fij = ', fij

        # define next term
        cfraction = (aux_sfSum - g[i])/vecpnorm(fi,p)
        # Recall that the coeff matrix satisfies vecpnorm(fi, p)=1 for most models.
        #cfraction = (aux_sfSum - g[i])

        if (cfraction > aux):
            aux = cfraction;

    #pX0 = 1 # asphericity is 'assumed' = 1
    # aux = aux*px0
    return aux



def H_distance_vertices_brute_force(X, Y, norm=2):
    r"""Compute Hausdorff distance between polytopes X and Y.

    The Hausdorff distance is defined as \max \{ \max_{x \in X}\min_{y \in Y} ||x-y||_p, \max_{y \in Y}\min_{x \in X} ||x-y||_p \} \}
    This algorithm is a straightforward implementation of the definition, that runs through all vertices.

    INPUTS:

    * "X, Y" : polyhedron objects.

    * "norm" : (default: 2). A value between 1 and infinity ('inf').

    OUTPUTS:

    * "dH" : Hausdorff distance between P and Q.

    """

    if X.base_ring() != Y.base_ring():
        raise ValueError('The base ring of the two input polytopes should be the same.')
    else:
        base_ring = X.base_ring()

    if norm == 'inf':
        norm = oo

    vertices_X = X.vertices_list()
    vertices_Y = Y.vertices_list()

    dXY = 0; dYX = 0

    # find dXY
    aux_ext = (vector(base_ring, vertices_X[0]) - vector(base_ring, vertices_Y[0])).norm(norm);
    for vx in vertices_X:
        aux_int = (vector(base_ring, vx) - vector(base_ring, vertices_Y[0])).norm(norm)
        for vy in vertices_Y:
            d = (vector(base_ring, vx) - vector(base_ring, vy)).norm(norm)
            if d < aux_int:
                aux_int = d

        if aux_int > aux_ext:
            aux_ext = aux_int
    dXY = aux_ext

    # find dYX
    aux_ext = (vector(base_ring, vertices_X[0]) - vector(base_ring, vertices_Y[0])).norm(norm);
    for vy in vertices_Y:
        aux_int = (vector(base_ring, vertices_X[0]) - vector(base_ring, vy)).norm(norm)
        for vx in vertices_X:
            d = (vector(base_ring, vx) - vector(base_ring, vy)).norm(norm)
            if d < aux_int:
                aux_int = d

        if aux_int > aux_ext:
            dYX = aux_int
    dYX = aux_ext

    return max(dXY, dYX)
