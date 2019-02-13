# Script Implementing the GGH Cryptosystem
# Reference: An Introduction to Mathematical Cryptography - Hoffstein, Pipher, Silverman

import numpy as np

def volume_lattice(a):
	'''
	Expects basis vectors as rows of matrix a
	'''
	return np.absolute(np.linalg.det(a))
def hadamard_ratio(a):
	'''
	Expects basis vectors as rows of matrix a
	'''
	vol = volume_lattice(a)
	dim = a.shape[0]
	norms_for_each_basis_vector = np.sum(np.abs(a)**2,axis=-1)**(1./2)
	product_of_norms = np.prod(norms_for_each_basis_vector)
	return (vol/product_of_norms)**(1./dim)

def gen_good_basis(d, dim):
	'''
	d: int param - coeffecients will vary between -d and d
	dim: dimension to be used for basis
	'''
	return np.random.randint(-d,d,(dim,dim))

# This function does not yield expected output
def gen_u():

	dim = 3
	initial_array = np.identity(dim)
	for i in xrange(0,10000):
		tmp = np.random.randint(0,1,(dim,dim))
		if np.any(tmp) > 0:
			initial_array = np.matmul(tmp,initial_array)
	print initial_array
# This function does not yield expected output
def gen_bad_basis(a):
	had_rat = 1
	dim = a.shape[0]
	initial_array = np.identity(dim)
	determinant = 1
	while (had_rat > 0.1) and (abs(determinant) != 1):
		tmp = np.random.randint(0,1,(dim,dim))
		initial_array = np.matmul(tmp,initial_array)
		had_rat = hadamard_ratio(initial_array)
		determinant = np.linalg.det(initial_array)
	print initial_array

# Encrypt a message given a public key, a message vector m and a random pertubation
def ggh_encrypt(public_mat,m, r):
	dim = public_mat.shape[0]
	#r = np.random.randint(-delta,delta,(dim,dim))
	# r commented for demo
	e = np.matmul(m,public_mat) + r
	return e

# Decrypt a message given a private key, an encrypted message and the associated public key
def ggh_decrypt(private_basis,e, public_basis):
	lam = np.matmul(e,np.linalg.inv(private_basis)) # express encrypted message as a combination of private basis
	lam = np.rint(lam) # round to nearest integer
	v = np.matmul(lam,private_basis) # Find a lattice vector close to e
	m = np.matmul(v,np.linalg.inv(public_basis)) #Express v as a linear combination of public basis, revealing m
	print np.rint(m)

def main():
	# Example from reference text
	priv_basis = np.array([[-97,19,19],[-36,30,86], [-184,-64,78]])
	pub_basis = np.array([[-4179163, -1882253,583183], [-3184353,-1434201, 444361],[-5277320,-2376852,736426]])
	m = np.array([86,-35,-32])
	r = np.array([-4,-3,2])
	ggh_decrypt(priv_basis,ggh_encrypt(pub_basis, m, r), pub_basis)
	#gen_u()

if __name__ == '__main__':
	main()
	