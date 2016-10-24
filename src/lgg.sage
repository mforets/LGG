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

    OUTPUTS:

    * "flowpipe"

    """

    # ---
    # Parse input
    # ---

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
    elif 'sage.geometry.polyhedron' not in str(type(X0)):
        # If X0 is not some type of polyhedron, try to set an initial point
        if type(X0) == list:
            X0 = Polyhedron(vertices = [X0], base_ring = base_ring)
        elif X0.is_vector():
            X0 = Polyhedron(vertices = [X0], base_ring = base_ring)
        else:
            raise ValueError('Initial state X0 not understood.')

    if B is None:
        # the system is homogeneous: dx/dt = Ax
        got_homogeneous = True
        raise NotImplementedError('Homogeneous system not implemented.')
    else:
        got_homogeneous = False
        if U is None:
            raise ValueError('Input range U is missing.')

    tau = kwargs['time_step'] if 'time_step' in kwargs else 1e-2

    t0 = kwargs['initial_time'] if 'initial_time' in kwargs else 0

    T = kwargs['time_horizon'] if 'time_horizon' in kwargs else 1

    N = kwargs['number_of_time_steps'] if 'number_of_time_steps' in kwargs else ceil(T/tau)

    directions = kwargs['directions'] if 'directions' in kwargs else {'select':'box'}

    # compute range of the input under B, V = BU
    V = polyhedron_linear_map(B, U, base_ring = base_ring)

    # ---
    # Compute first element of the approximating sequence
    # ---

    # compute matrix exponential exp(A*tau)
    Phi_tau = expm(np.multiply(A, tau))

    # compute exp(tau*A)X0
    expX0 = polyhedron_linear_map(Phi_tau, X0, base_ring = base_ring)

    # compute the initial over-approximation
    tau_V = polyhedron_linear_map(tau*np.identity(n), V)

    # compute the bloating factor
    Ainfty = matrix_sup_norm(A)
    RX0 = polyhedron_sup_norm(X0)
    RV = polyhedron_sup_norm(V)

    unitBall = BoxInfty(center = zero_vector(n), radius = 1, base_ring = base_ring)
    alpha_tau = (exp(tau*Ainfty) - 1 - tau*Ainfty)*(RX0 + RV/Ainfty)
    alpha_tau_B = polyhedron_linear_map(alpha_tau*np.identity(n), unitBall)

    # compute the first element of the approximating sequence, Omega_0
    aux = expX0.Minkowski_sum(tau_V)
    Omega0 = X0.convex_hull(aux.Minkowski_sum(alpha_tau_B))

    # ---
    # Build the sequence of approximations Omega_i
    # ---

    # build the Omega_i
    beta_tau = (exp(tau*Ainfty) - 1 - tau*Ainfty)*(RV/Ainfty)
    beta_tau_B = polyhedron_linear_map(beta_tau*np.identity(n), unitBall)
    W_tau = tau_V.Minkowski_sum(beta_tau_B)

    # receives a direction d, and outputs the support function rho_i of Omega_i at d
    def _Omega_i_supports(d):
        r = []
        s = []
        rhoi = []

        r.append(d)
        s.append(0)
        rhoi.append(supp_fun_polyhedron(Omega0, d))

        for i in [0..N-2]:
            r.append(np.dot(Phi_tau.transpose(),r[i]))
            s.append(s[i] + supp_fun_polyhedron(W_tau, r[i]))
            rhoi.append(s[i+1] + supp_fun_polyhedron(Omega0, r[i+1]))

        return rhoi

    # ---
    # Define the directions
    # ---

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


    # directions as numpy array, and number of directions
    dArray = np.array(dList)
    k = len(dArray)

    # compute family of support functions
    Omega_i_Family_SF = list()
    for i in range(len(dArray)):
        d = dArray[i]
        Omega_i_Family_SF.append( _Omega_i_supports(d) )


    # ---
    # Compute the flowpipe, all the elements of the approximating sequence
    # ---

    # build each polyhedron using the support functions over-approximation
    # we have N polyhedrons
    Omega_i_Poly = list()

    #run over polytopes
    for i in range(N):

        # for each polyhedra, I have to use all directions that I know
        A = matrix(RR,k,n); b = vector(RR,k)
        for j in range(k): #run over directions
            s_fun = Omega_i_Family_SF[j][i]
            A.set_row(j, dList[j])
            b[j] = s_fun

        Omega_i_Poly.append( PolyhedronFromHSpaceRep(A, b, base_ring = base_ring) )

    # return Omega0 ?
    # return [Omega0] + Omega_i_Poly

    return Omega_i_Poly



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
