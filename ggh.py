# Script Implementing the GGH Cryptosystem
# Reference: An Introduction to Mathematical Cryptography - Hoffstein, Pipher, Silverman

import sympy
from lattice import *
import matplotlib.pyplot as plt



def gen_good_basis(n, scale, *, ratio=None):
    """
    Generate a "good" basis for a random lattice with a large Hadamard ratio.
    
    Parameters:
        n - the number of dimensions of the good basis
        scale - range of random values of the good basis
        ratio - (default: 0.95) the minimum hadamard ratio of the good basis
    """
    
    if ratio is None:
        ratio = 0.95
    
    basis = Basis(random_matrix(n, scale=scale))
    while hadamard_ratio(basis) < ratio:
        basis = gram_schmidt(Basis(random_matrix(n, scale=scale))).map(round)
    
    return basis



def gen_bad_basis(basis, *, ratio=None, min_ratio=None, scale=None):
    """
    Generate a "bad" basis with a small Hdamard ratio and large entries that
    spans the same lattice as a "good" input lattice.
    
    Parameters:
        basis - the original basis for a lattice
        ratio - (default: 0.05) the maximum hadamard ratio of the bad basis
        scale - (default: 1000) the approximate scale factor of the bad basis
    """
    
    if ratio is None:
        ratio = 0.05
    
    if min_ratio is None:
        min_ratio = 0.00
    
    if scale is None:
        scale = 100
    
    n = len(basis)
    
    new_basis = Basis(basis.mat() * random_unimodular_matrix(n, scale=scale))
    while hadamard_ratio(new_basis) > ratio or hadamard_ratio(new_basis) < min_ratio:
        new_basis = Basis(basis.mat() * random_unimodular_matrix(n, scale=scale))
    
    return new_basis



def ggh_encrypt(public_basis, message, *, delta=None, rand=None):
    """
    Returns an encrypted message using the GGH encryption algorithm.
    
    Parameters:
        public_basis - the public key used in GGH (a bad basis for a lattice)
        message - a vector containing the unencrypted message
        delta - (default: 10) the range of values for the random vector perturbation
        rand - (randomized by default) the random vector perturbation
    """
    
    if delta is None:
        delta = 10
    
    if rand is None:
        rand = random_matrix(len(public_basis), 1, scale=delta)
    
    return public_basis.mat() * message + rand



def ggh_decrypt(private_basis, public_basis, encrypted_message):
    """
    Returns a decrypted message using the GGH encryption algorithm.
    
    Parameters:
        private_basis - the private key used in GGH (a good basis for a lattice)
        public_basis  - the public key used in GGH (a bad basis for a lattice)
        message - a vector containing the encrypted message
    """
    
    # Find the lattice vector v closest to the vector of the encrypted message
    v = private_basis.mat() * (private_basis.inv() * encrypted_message).applyfunc(round)
    
    # Express v as a linear combination of public basis, revealing the original message
    return (public_basis.inv() * v).applyfunc(round)



def ggh_crack(public_basis, encrypted_message):
    """
    Returns a decrypted message by attempting to generate a good basis using the
    LLL algorithm. If the basis returned by LLL is not sufficiently small or
    orthogonal, then the decrypted message is likely to be incorrect.
    """
    
    # Find a good basis for the lattice
    good_basis = lll_lattice_reduction(public_basis)
    
    # Decrypt using the good basis as if it's the private basis
    return ggh_decrypt(good_basis, public_basis, encrypted_message)



def run_example_ggh():
    """
    Runs the example GGH from reference text
    """
    
    private_basis = Basis([
        [-97, 19, 19],
        [-36, 30, 86],
        [-184, -64, 78]
    ], row_vectors=True)
    
    public_basis = Basis([
        [-4179163, -1882253, 583183],
        [-3184353, -1434201, 444361],
        [-5277320, -2376852, 736426]
    ], row_vectors=True)
    
    rand = Matrix([-4, -3, 2])
    message = Matrix([86, -35, -32])
    print('Original message:')
    print(message)
    
    encrypted = ggh_encrypt(public_basis, message, rand=rand)
    print()
    print('Encrypted message:')
    print(encrypted)
    
    decrypted = ggh_decrypt(private_basis, public_basis, encrypted)
    print()
    print('Decrypted message:')
    print(decrypted)



def run_random_ggh():
    private_basis = gen_good_basis(3, 100)
    public_basis = gen_bad_basis(private_basis)
    
    print()
    print(private_basis)
    print(hadamard_ratio(private_basis))
    
    print()
    print(public_basis)
    print(hadamard_ratio(public_basis))



def run_crack_ggh():
    public_basis = Basis([
        [-4179163, -1882253, 583183],
        [-3184353, -1434201, 444361],
        [-5277320, -2376852, 736426],
    ], row_vectors=True)
    
    original = Matrix([86, -35, -32])
    encrypted = Matrix([-79081427, -35617462, 11035473])
    decrypted = ggh_crack(public_basis, encrypted)
    
    print('decrypted:', decrypted)
    print('original:', original)



def smooth_data(x, y, factor=5):
    n = len(x)
    m = n
    d = (x[-1] - x[0]) / (m - 1)
    
    def w(i, j):
        return max(0, factor * d - abs(xp[i] - x[j]))
    
    xp = [0 for i in range(m)]
    yp = [0 for i in range(m)]
    
    for i in range(m):
        xp[i] = x[i]
        yp[i] = sum(y[j] * w(i, j) for j in range(n)) / sum(w(i, j) for j in range(n))
    
    return xp, yp



def test_ggh_decryption():
    """
    Test how often GGH decryption fails for different Hadamard ratios of the
    private basis. How does this vary with n?
    """
    
    scale = 100
    trials = 128
    
    data = []
    for n in range(4, 8):
        
        points = []
        for i in range(trials):
            private = gen_good_basis(n, scale, ratio=0.5)
            public = gen_bad_basis(private, ratio=0.05)
            ratio = hadamard_ratio(private)
            
            msg = random_matrix(n, 1, scale=scale)
            rec = ggh_decrypt(private, public, ggh_encrypt(public, msg, delta=10))
            
            if rec == msg:
                points.append((ratio, 1))
            else:
                points.append((ratio, 0))
        
        points.sort(key=lambda x: x[0])
        print(i, points)
        
        ratios = [x[0] for x in points]
        success = [x[1] for x in points]
        data.append((n, ratios, success))
    
        #for n, ratios, success in data:
        plt.clf()
        fig = plt.figure()
        
        plt.title('Decrypting GGH Success Rate (n={})'.format(n))
        plt.xlabel('Hadamard ratio of the private basis')
        plt.ylabel('Rate of success (decrypting an encrypted message)')
        plt.xlim(0.5, 1.0)
        plt.ylim(-0.2, 1.2)
        plt.plot(*smooth_data(ratios, success, 16))
        
        fig.savefig('../decrypting_success_{}.png'.format(n), dpi=fig.dpi)
    
    return data



def test_ggh_cracking():
    """
    Test how often GGH cracking fails for different Hadamard ratios of the
    public basis. How does this vary with n?
    """
    
    scale = 100
    trials = 128
    
    data = []
    for n in range(4, 8):
        
        points = []
        for i in range(trials):
            private = gen_good_basis(n, scale, ratio=0.8)
            public = gen_bad_basis(private, ratio=0.5)
            ratio = hadamard_ratio(public)
            
            msg = random_matrix(n, 1, scale=scale)
            rec = ggh_crack(public, ggh_encrypt(public, msg, delta=10))
            
            if rec == msg:
                points.append((ratio, 1))
            else:
                points.append((ratio, 0))
        
        points.sort(key=lambda x: x[0])
        print(i, points)
        
        ratios = [x[0] for x in points]
        success = [x[1] for x in points]
        data.append((n, ratios, success))
    
        #for n, ratios, success in data:
        plt.clf()
        fig = plt.figure()
        
        plt.title('Cracking GGH Success Rate (n={})'.format(n))
        plt.xlabel('Hadamard ratio of the public basis')
        plt.ylabel('Rate of success (cracking an encrypted message)')
        plt.xlim(0.0, 0.1)
        plt.ylim(-0.2, 1.2)
        plt.plot(*smooth_data(ratios, success, 16))
        
        fig.savefig('../cracking_success_{}.png'.format(n), dpi=fig.dpi)
    
    return data



def main():
    print(test_ggh_decryption())
    #run_crack_ggh()
    #run_example_ggh()
    #run_random_ggh()



if __name__ == '__main__':
    main()

