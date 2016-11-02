r"""

TO-DO:

* Check if it is relevant to Cythonize some parts. What can be done with a user-defined
class (not native C). For instance, with the Polyhedron?

* Use an internal supp_fun_polyhedron, because W_tau and Omega_0 do not change, so we may avoid some overhead
due to function calls and preparation of the LP. In experiments, Gurobi is being slower. For instance, I got 8s in GLPK
against 14s in Gurobi. I think this may be from unnecessary overhead.

* Understand why in some configurations it does not work. For instance, if X0 is QQ and base_ring is set to QQ. However,
if X0 is RDF and base_ring is set to QQ is does work (only with GLPK). Also, it would be intersting to make it work with Gurobi.
Currently, if X0 is RDF and base_ring is set to QQ, then Gurobi will give infeas or unbdd problem.

* Homogeneous case is still not correct. To see this, go to
simple2d.ipynb
We expect to get the same results as with mu -> 0.

"""
import sys
sys.path.append('..')

from lib.polyFunctions_core import *
from lib.norms import matrix_sup_norm

import numpy as np
from scipy.linalg import expm, sinm, cosm
import random

class NotImplementedException(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

def compute_flowpipe(A=None, X0=None, B=None, U=None, **kwargs):
    r"""Implements LGG reachability algorithm for the linear continuous system dx/dx = Ax + Bu.

    INPUTS:

    * "A" - The coefficient matrix of the system.

    * "X0" - The initial state.

    * "B" -

    * "U" -

    * "time_step" - (default = 1e-2). Time step.

    * "initial_time" - (default = 0). The initial time.

    * "time_horizon" - (default = 1). The final time.

    * "number_of_time_steps" - (default = ceil(T/tau)). Number of time steps.

    * "directions" - (default: random, and a box). A dictionary.

    * "solver" - LP solver. Valid options are:
        * 'GLPK' (default).
        * 'Gurobi' - not tested.

    * "base_ring" - Base ring where polyhedral computations are performed.
        Valid options are:
        * QQ - (default) rational field
        * RDF - real double field

    OUTPUTS:

    * "flowpipe"

    """

    # ################
    # Parse input    #
    # ################
    if A is None:
        raise ValueError('System matrix A is missing.')
    else:
        if 'sage.matrix' in str(type(A)):
            n = A.ncols()
        elif type(A) == np.ndarray:
            n = A.shape[0]

    base_ring = kwargs['base_ring'] if 'base_ring' in kwargs else QQ

    if X0 is None:
        raise ValueError('Initial state X0 is missing.')
    elif 'sage.geometry.polyhedron' not in str(type(X0)) and type(X0) == list:
        # If X0 is not some type of polyhedron, set an initial point
        X0 = Polyhedron(vertices = [X0], base_ring = base_ring)
    elif 'sage.geometry.polyhedron' not in str(type(X0)) and X0.is_vector():
        X0 = Polyhedron(vertices = [X0], base_ring = base_ring)
    elif 'sage.geometry.polyhedron' in str(type(X0)):
        # ensure that all input sets are on the same ring
        # not sure about this
        if 1==0:
            if X0.base_ring() != base_ring:
                [F, g] = PolyhedronToHSpaceRep(X0)
                X0 = PolyhedronFromHSpaceRep(F, g, base_ring=base_ring)
    else:
        raise ValueError('Initial state X0 not understood')

    if B is None:
        # the system is homogeneous: dx/dt = Ax
        got_homogeneous = True
    else:
        got_homogeneous = False
        if U is None:
            raise ValueError('Input range U is missing.')

    tau = kwargs['time_step'] if 'time_step' in kwargs else 1e-2

    t0 = kwargs['initial_time'] if 'initial_time' in kwargs else 0

    T = kwargs['time_horizon'] if 'time_horizon' in kwargs else 1

    global N
    N = kwargs['number_of_time_steps'] if 'number_of_time_steps' in kwargs else ceil(T/tau)

    directions = kwargs['directions'] if 'directions' in kwargs else {'select':'box'}

    global solver
    solver = kwargs['solver'] if 'solver' in kwargs else 'GLPK'

    global verbose
    verbose = kwargs['verbose'] if 'verbose' in kwargs else 0

    # #######################################################
    # Compute first element of the approximating sequence   #
    # #######################################################

    global Phi_tau, Omega0, W_tau

    if got_homogeneous: # dx/dx = Ax

        # compute range of the input under B, V = BU
        #V = polyhedron_linear_map(B, U, base_ring = base_ring)

        # compute matrix exponential exp(A*tau)
        Phi_tau = expm(np.multiply(A, tau))

        # compute exp(tau*A)X0
        expX0 = polyhedron_linear_map(Phi_tau, X0, base_ring = base_ring)

        # compute the initial over-approximation
        #tau_V = polyhedron_linear_map(tau*np.identity(n), V)

        # compute the bloating factor
        Ainfty = matrix_sup_norm(A)
        RX0 = polyhedron_sup_norm(X0)
        #RV = polyhedron_sup_norm(V)

        unitBall = BoxInfty(center = zero_vector(n), radius = 1, base_ring = base_ring)
        alpha_tau = (exp(tau*Ainfty) - 1 - tau*Ainfty)*(RX0)
        alpha_tau_B = polyhedron_linear_map(alpha_tau*np.identity(n), unitBall, base_ring = base_ring)

        # compute the first element of the approximating sequence, Omega_0
        #aux = expX0.Minkowski_sum(tau_V)
        Omega0 = X0.convex_hull(alpha_tau_B)

        #beta_tau = (exp(tau*Ainfty) - 1 - tau*Ainfty)*(RV/Ainfty)
        #beta_tau_B = polyhedron_linear_map(beta_tau*np.identity(n), unitBall)
        W_tau = Polyhedron(vertices = [], ambient_dim=n) # NOT TESTED
        # if W_tau = [], then supp_fun_polyhedron is set to return 0
        #W_tau = tau_V.Minkowski_sum(beta_tau_B)

    else: # dx/dx = Ax + Bu

        # compute range of the input under B, V = BU
        V = polyhedron_linear_map(B, U, base_ring = base_ring)

        # compute matrix exponential exp(A*tau)
        Phi_tau = expm(np.multiply(A, tau))

        # compute exp(tau*A)X0
        expX0 = polyhedron_linear_map(Phi_tau, X0, base_ring = base_ring)

        # compute the initial over-approximation
        tau_V = polyhedron_linear_map(tau*np.identity(n), V, base_ring = base_ring)

        # compute the bloating factor
        Ainfty = matrix_sup_norm(A)
        RX0 = polyhedron_sup_norm(X0)
        RV = polyhedron_sup_norm(V)

        unitBall = BoxInfty(center = zero_vector(n), radius = 1, base_ring = base_ring)
        alpha_tau = (exp(tau*Ainfty) - 1 - tau*Ainfty)*(RX0 + RV/Ainfty)
        alpha_tau_B = polyhedron_linear_map(alpha_tau*np.identity(n), unitBall, base_ring = base_ring)

        # compute the first element of the approximating sequence, Omega_0
        aux = expX0.Minkowski_sum(tau_V)
        Omega0 = X0.convex_hull(aux.Minkowski_sum(alpha_tau_B))

        beta_tau = (exp(tau*Ainfty) - 1 - tau*Ainfty)*(RV/Ainfty)
        beta_tau_B = polyhedron_linear_map(beta_tau*np.identity(n), unitBall, base_ring = base_ring)
        W_tau = tau_V.Minkowski_sum(beta_tau_B)

    # #######################################################
    # Generate directions                                   #
    # #######################################################
    if directions['select'] == 'box':

        if n==2:
            theta = [0,pi/2,pi,3*pi/2] # box
            dList = [vector(RR,[cos(t), sin(t)]) for t in theta]

        else: # directions of hypercube
            dList = []
            dList += [-identity_matrix(n).column(i) for i in range(n)]
            dList += [identity_matrix(n).column(i) for i in range(n)]

    elif directions['select'] == 'oct':

        if n != 2:
            raise NotImplementedError('Directions select octagon not implemented for n other than 2. Try box.')

        theta = [i*pi/4 for i in [0..7]] # octagon
        dList = [vector(RR,[cos(t), sin(t)]) for t in theta]

    elif directions['select'] == 'random':

        order = directions['order'] if 'order' in directions else 12

        if n == 2:
            theta = [random.uniform(0, 2*pi.n(digits=5)) for i in range(order)]
            dList = [vector(RR,[cos(theta[i]), sin(theta[i])]) for i in range(order)]
        else:
            raise NotImplementedError('Directions select random not implemented for n greater than 2. Try box.')

    elif directions['select'] == 'custom':

        dList = directions['dList']

    else:

        raise TypeError('Directions dictionary not understood.')


    # transform directions to numpy array, and get number of directions
    dArray = np.array(dList)
    k = len(dArray)

    # compute family of support functions
    Omega_i_Family_SF = list()
    for i in range(len(dArray)):
        d = dArray[i]
        Omega_i_Family_SF.append( _Omega_i_supports(d) )

    # ################################################
    # Build the sequence of approximations Omega_i   #
    # ################################################

    # each polyhedron is built using the support functions over-approximation
    # we have N polyhedrons
    Omega_i_Poly = list()

    # This loop can be parallelized
    for i in range(N):    # run over polyhedra

        # for each polyhedron, use all directions
        A = matrix(RR,k,n); b = vector(RR,k)

        for j in range(k): #run over directions
            s_fun = Omega_i_Family_SF[j][i]
            A.set_row(j, dList[j])
            b[j] = s_fun

        Omega_i_Poly.append( PolyhedronFromHSpaceRep(A, b, base_ring = base_ring) )

    # return Omega0 ?
    # return [Omega0] + Omega_i_Poly

    return Omega_i_Poly

def _Omega_i_supports(d):
    r"""Receives a direction d, and outputs the support function
    rho_i of Omega_i at d"""

    r = []
    s = []
    rhoi = []

    r.append(d)
    s.append(0)
    rhoi.append(supp_fun_polyhedron(Omega0, d, solver=solver, verbose=verbose))

    for i in [0..N-2]:
        r.append(np.dot(Phi_tau.transpose(),r[i]))
        #print r[i]
        #print W_tau.base_ring()
        #print W_tau.inequalities_list()
        #print supp_fun_polyhedron(W_tau, r[i], solver=solver, verbose=verbose)
        #print '===================='
        s.append(s[i] + supp_fun_polyhedron(W_tau, r[i], solver=solver, verbose=verbose))
        rhoi.append(s[i+1] + supp_fun_polyhedron(Omega0, r[i+1], solver=solver, verbose=verbose))

    return rhoi

def plot_flowpipe(fp, **kwargs):

    from sage.geometry.polyhedron.plot import Projection

    n = fp[0].ambient_dim()

    if n == 2:

        # plot the result
        myFig = Graphics()
        myFig = sum(p.plot(alpha=0.5) for p in fp)
        myFig.show()

    elif n>2:

        v = kwargs['projection_directions']

        fp_proj = [Projection(p, proj = lambda x : [ x[v[0]], x[v[1]] ]) for p in fp]

        # plot the result
        myFig = Graphics()
        myFig = sum(p.plot() for p in fp_proj)
        myFig.show()

    return
