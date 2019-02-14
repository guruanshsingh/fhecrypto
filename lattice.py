
"""
This module provides useful methods for working with lattices, especially for
cryptographic purposes. Many of the algorithms included here are taken from
"An Introduction to Mathematical Cryptography" (the textbook) by Hoffstein,
Pipher, and Silverman. 
"""



import random

import sympy
from sympy import Matrix



def prod(iterable, *, empty=1):
    """
    Returns the left-associative product of the elements in an iterable.
    
    Paramters:
        iterable - an iterable whose elements can be multiplied together
        empty    - (default: 1) the value returned for an empty product
    """
    
    p = empty
    for x in iterator:
        p *= x
    return p



def rows(matrix):
    """
    Returns a list of the rows of a matrix in order from top to bottom.
    """
    
    return [matrix[i,:] for i in range(matrix.shape[0])]



def columns(matrix):
    """
    Returns a list of the columns of a matrix in order from left to right.
    """
    
    return [matrix[:,i] for i in range(matrix.shape[1])]



def lattice_volume(basis):
    """
    Returns the volume of the fundamental domain of a lattice.
    
    Parameters:
        basis - a list of vectors that span a lattice
    """
    
    return abs(Matrix.hstack(*basis).det())



def random_unimodular_matrix(n, *, scale=1000):
    """
    Returns a random n by n unimodular matrix. This is done by multiplying
    together random upper and lower triangular matrices whose diagonals contain
    only 1 and -1.
    
    Parameters:
        n     - the size of the resulting matrix
        scale - (defualt: 1000) the scale of entries in the resulting matrix
    """
    
    m = Matrix.eye(n)
    r = Matrix.eye(n)
    
    while max(m) <= scale:
        # Randomize the upper triangle of matrix r
        for i in range(n):
            r[i,i] = random.choice([-1, 1])
            for j in range(i + 1, n):
                r[i,j] = random.randint(-1, 1)
        
        # Randomly alternate between upper and lower triangular matrices
        if random.choice([False, True]):
            m *= r
        else:
            m *= r.T
    
    return m



def hadamard_ratio(basis):
    """
    Returns the Hadamard ratio of a basis of vectors
    
    Parameters:
        basis - a list of vectors that span a lattice
    """
    
    n = len(basis)
    vol = lattice_volume(basis)
    product_of_norms = prod(v.norm() for v in basis)
    return (vol / product_of_norms) ** (1 / n)



def gram_schmidt(basis):
    """
    Returns a basis of pairwise orthogonal vectors which span the same space as
    the original basis.
    
    Parameters:
        basis - a list of linearly independent vectors
    """
    
    new_basis = []
    for v in basis:
        for u in new_basis:
            v = v - (u.dot(v) / u.norm() ** 2) * u
        new_basis.append(v)
    
    return new_basis



def gaussian_lattice_reduction(u, v):
    """
    Returns the two smallest vectors which span a 2-dimensional lattice.
    
    This algorithm is described in section 7.13.1 (page 436) of the textbook.
    
    Parameters:
        u - the first basis vector of a 2D lattice
        v - the second basis vector of a 2D lattice
    """
    
    while True:
        if u.norm() ** 2 > v.norm() ** 2:
            u, v = v, u
        
        m = round(u.dot(v) / u.norm() ** 2)
        
        if m == 0:
            return [u, v]
        
        v = v - m * u



def lll_lattice_reduction(basis):
    """
    Returns a list of "small" vectors which span the same lattice as the
    original basis.
    
    This algorithm is described in section 7.13.2 (page 439) of the textbook.
    
    Parameters:
        basis - a list of linearly independent vectors
    """
    
    basis = basis[:]
    n = len(basis)
    
    def mu(i, j):
        return basis[i].dot(tmp[j]) / tmp[j].norm() ** 2
    
    k = 1
    while k < n:
        tmp = gram_schmidt(basis)
        
        for j in range(k - 1, -1, -1):
            basis[k] = basis[k] - basis[j] * round(mu(k, j))
        
        if 4 * (tmp[k].norm() / tmp[k - 1].norm()) ** 2 >= 3 - 4 * mu(k, k - 1) ** 2:
            k = k + 1
        else:
            basis[k], basis[k - 1] = basis[k - 1], basis[k]
            k = max(k - 1, 1)
    
    return basis



