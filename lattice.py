
"""
This module provides useful methods for working with lattices, especially for
cryptographic purposes. Many of the algorithms included here are taken from
"An Introduction to Mathematical Cryptography" (the textbook) by Hoffstein,
Pipher, and Silverman. 
"""



import random

import sympy
from sympy import Matrix



class Basis:
    """
    Represents the basis of a vector-space or lattice.
    """
    
    def __init__(self, arg):
        if isinstance(arg, Basis):
            self.vectors = [vec.copy() for vec in arg]
            
        elif isinstance(arg, sympy.MatrixBase):
            self.vectors = [arg[:, i] for i in range(arg.cols)]
            
        elif isinstance(arg, list):
            self.vectors = [Matrix(x) for x in arg]
            
        else:
            try:
                self.vectors = [Matrix(x) for x in arg]
            except TypeError:
                raise TypeError('Cannot construct a Basis from an object of type {}'.format(type(arg)))
    
    def __str__(self):
        return self.mat().__str__()
    
    def __repr__(self):
        return 'Basis({})'.format(self.vectors.__str__())
    
    def __len__(self):
        return self.vectors.__len__()
    
    def __getitem__(self, i):
        return self.vectors.__getitem__(i)
    
    def __setitem__(self, i, val):
        return self.vectors.__setitem__(i, val)
    
    def __iter__(self):
        return self.vectors.__iter__()
    
    def copy(self):
        return Basis([vec.copy() for vec in self.vectors])
    
    def vol(self):
        return abs(self.mat().det())
    
    def mat(self):
        return Matrix.hstack(*self.vectors)
    
    def inv(self):
        return self.mat().inv()



def prod(iterable, *, empty=1):
    """
    Returns the left-associative product of the elements in an iterable.
    
    Paramters:
        iterable - an iterable whose elements can be multiplied together
        empty    - (default: 1) the value returned for an empty product
    """
    
    p = empty
    for x in iterable:
        p *= x
    return p



def hadamard_ratio(basis):
    """
    Returns the Hadamard ratio of a basis of vectors
    
    Parameters:
        basis - a list of vectors that span a lattice
    """
    
    product_of_norms = prod(v.norm() for v in basis)
    if product_of_norms == 0:
        return 0
    
    return (basis.vol() / product_of_norms) ** (1 / len(basis))



def random_matrix(n, m=None, *, scale=None):
    """
    Returns a random n by n unimodular matrix. The matrix should be composed of
    elements on the order of 'scale'.
    
    Parameterss:
        n     - the number of rows of the matrix
        m     - (default: n) the number of columns of the matrix
        scale - (default: 1000) the scale of entries in the matrix
    """
    
    if m is None:
        m = n
    
    if scale is None:
        scale = 1000
    
    return Matrix([[random.randint(-scale, scale) for c in range(m)] for r in range(n)])

def random_unimodular_matrix(n, *, scale=None, rand=None):
    """
    Returns a random n by n unimodular matrix. The matrix should be composed of
    elements on the order 'scale'. The matrix is calculared by multiplying
    together random upper and lower triangular matrices.
    
    Parameters:
        n     - the size of the matrix
        scale - (defualt: 1000) the scale of entries in the matrix
        rand  - (default: (-1, 1)) the range of triangular matrix entries
    """
    
    if scale is None:
        scale = 1000
    
    if rand is None:
        rand = (-1, 1)
    
    m = Matrix.eye(n)
    r = Matrix.eye(n)
    
    while max(m) <= scale:
        # Randomly permute the columns of m
        for i in range(n):
            m.col_swap(i, random.randrange(i, n))
        
        # Randomize the upper triangle of matrix r
        for i in range(n):
            for j in range(i + 1, n):
                r[i,j] = random.randint(*rand)
        
        # Randomly alternate between upper and lower triangular matrices
        if random.choice([False, True]):
            m *= r
        else:
            m *= r.T
    
    return m



def projection_factor(u, v):
    """
    Returns the scaling factor of vector u projected onto vector v.
    """
    
    return u.dot(v) / v.dot(v)



def projection(u, v):
    """
    Returns the projection of vector u onto vector v.
    """
    
    return (u.dot(v) / v.dot(v)) * v



def gram_schmidt(basis):
    """
    Returns a basis of pairwise orthogonal vectors which span the same space as
    the original basis.
    
    Parameters:
        basis - a list of linearly independent vectors
    """
    
    new_basis = basis.copy()
    for i in range(len(basis)):
        for j in range(i):
            new_basis[i] -= projection(new_basis[i], basis[j])
    
    return new_basis



def gaussian_lattice_reduction(basis):
    """
    Returns the two smallest vectors which span a 2-dimensional lattice.
    
    This algorithm is described in section 7.13.1 (page 436) of the textbook.
    
    Parameters:
        u - the first basis vector of a 2D lattice
        v - the second basis vector of a 2D lattice
    """
    
    if len(basis) != 2:
        raise ValueError('Gaussian reduction can only be done on 2-dimensional lattices')
    
    u, v = basis
    
    while True:
        if u.norm() > v.norm():
            u, v = v, u
        
        m = round(projection_factor(u, v))
        
        if m == 0:
            return Basis([u, v])
        
        v = v - m * u



def lll_lattice_reduction(basis):
    """
    Returns a list of "small" vectors which span the same lattice as the
    original basis.
    
    This algorithm is described in section 7.13.2 (page 439) of the textbook.
    
    Parameters:
        basis - a list of linearly independent vectors
    """
    
    n = len(basis)
    basis = basis.copy()
    
    def mu(i, j):
        return projection_factor(basis[i], ortho_basis[j])
    
    k = 1
    while k < n:
        ortho_basis = gram_schmidt(basis)
        
        for j in range(k - 1, -1, -1):
            basis[k] = basis[k] - basis[j] * round(mu(k, j))
        
        if 4 * (ortho_basis[k].norm() / ortho_basis[k - 1].norm()) ** 2 >= 3 - 4 * mu(k, k - 1) ** 2:
            k = k + 1
        else:
            basis[k], basis[k - 1] = basis[k - 1], basis[k]
            k = max(k - 1, 1)
    
    return basis



