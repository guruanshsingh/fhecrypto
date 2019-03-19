
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
    
    def __init__(self, arg, *, row_vectors=False):
        try:
            if isinstance(arg, Basis):
                self.matrix = arg.matrix.copy()
            elif isinstance(arg, list) and isinstance(arg[0], Matrix):
                self.matrix = Matrix.hstack(*arg)
            elif isinstance(arg, tuple) and isinstance(arg[0], Matrix):
                self.matrix = Matrix.hstack(*arg)
            else:
                self.matrix = Matrix(arg)
            
            if row_vectors:
                self.matrix = self.matrix.T
            
        except TypeError:
            raise TypeError('Cannot construct a Basis from an object of type {}'.format(type(arg)))
    
    def __str__(self):
        return self.matrix.__str__()
    
    def __repr__(self):
        return 'Basis({})'.format(self.matrix.__str__())
    
    def __len__(self):
        return self.matrix.cols
    
    def __getitem__(self, i):
        return self.matrix.__getitem__((slice(None), i))
    
    def __setitem__(self, i, val):
        return self.matrix.__setitem__((slice(None), i), val)
    
    def __iter__(self):
        for i in range(self.__len__()):
            yield self.__getitem__(i)
    
    def copy(self):
        return Basis(self.matrix.copy())
    
    def vol(self):
        return abs(self.matrix.det())
    
    def mat(self):
        return self.matrix
    
    def inv(self):
        return self.matrix.inv()
    
    def map(self, f):
        return Basis(self.matrix.applyfunc(f))



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
            new_basis[i] -= projection(basis[i], new_basis[j])
    
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



def lll_lattice_reduction(basis, lovasz_const=None):
    """
    Returns a list of "small" vectors which span the same lattice as the
    original basis.
    
    This algorithm is described in section 7.13.2 (page 439) of the textbook.
    
    Parameters:
        basis - a list of linearly independent vectors
        lovasz_const - (default: 3/4) the constant used by the lovasz condition
    """
    
    if lovasz_const is None:
        lovasz_const = 3/4
    
    n = len(basis)
    basis = basis.copy()
    
    def mu(i, j):
        return projection_factor(basis[i], ortho_basis[j])
    
    k = 1
    while k < n:
        ortho_basis = gram_schmidt(basis)
        
        for j in reversed(range(0, k)):
            basis[k] = basis[k] - basis[j] * round(mu(k, j))
        
        if ortho_basis[k].norm() ** 2 >= (lovasz_const - mu(k, k - 1) ** 2) * ortho_basis[k - 1].norm() ** 2:
            k = k + 1
        else:
            basis[k], basis[k - 1] = basis[k - 1], basis[k]
            k = max(k - 1, 1)
    
    return basis



def main():
    b = Basis([
        [19,  2, 32, 46,  3, 33],
        [15, 42, 11,  0,  3, 24],
        [43, 15,  0, 24,  4, 16],
        [20, 44, 44,  0, 18, 15],
        [ 0, 48, 35, 16, 31, 31],
        [48, 33, 32,  9,  1, 29]
    ], row_vectors=True)
    
    print('Original basis:', b)
    print('LLL basis:', lll_lattice_reduction(b))



if __name__ == '__main__':
    main()

