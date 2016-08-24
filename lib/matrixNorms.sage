
#===================================================================================#===================================================================================
#Functions to compute p-norms

# (note that numpy has its own functions to compute norms of vectors and matrices) 

# Computing usual matrix norms (p=1,2,infty)
# MISSING: estimate on the matrix norm for arbitrary p (cf. N. Higham algorithms) 

#We compute the matrix norm induced by vector p-norms, only for p=1,2,infty. The general case is left for future work. Hints are found in:
#http://mathoverflow.net/questions/39148/efficiently-computing-a-matrixs-induced-p-norm

#We compute the vector p norm for 1 <=p < infty with vecpnorm(b,p) and vecinftynorm(b) for p = infty
#===================================================================================#===================================================================================
import numpy as np

def mat1norm(A):
    aux = 0
    for j in range(A.ncols()):
        aux2=0
        for i in range(A.nrows()):
            aux2 += abs(A[i,j])
        if (aux2 > aux):
            aux = aux2

    return aux

def matinftynorm(A):

    try:
        [nrows, ncols] = A.shape
    except:
        nrows = A.nrows(); ncols = A.ncols()

    aux = 0
    for i in range(nrows):
        aux2=0
        for j in range(ncols):
            aux2 += abs(A[i,j])
        if (aux2 > aux):
            aux = aux2

    return aux

def mat2norm(A):
    Q = A.H*A
    eval, evec = np.linalg.eig(Q.change_ring(RR))
    aux = sqrt(np.max(eval))

    return aux


def vec1norm(b):
    aux = 0
    for i in range(len(b)):
        aux += abs(b[i])

    return aux


def vecinftynorm(b):
    aux = 0
    for i in range(len(b)):
        if (abs(b[i]) > aux):
            aux = abs(b[i]) 

    return aux


def vec2norm(b):
    aux = 0
    for i in range(len(b)):
        aux += b[i]^2

    return sqrt(aux)



# p is a number: 1 <= p < infty
def vecpnorm(b,p): 
    aux = 0
    for i in range(len(b)):
        aux += abs(b[i])^p

    return aux^(1/p)









