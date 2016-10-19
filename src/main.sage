import sys
sys.path.insert(1, '..')

from lib.polyFunctions_core import PolyhedronToHSpaceRep, PolyhedronFromHSpaceRep

def computePoly(A, b):

    P = PolyhedronFromHSpaceRep(A, b, base_ring = QQ)
    
    return P