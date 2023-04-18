import re
import time
import numpy as np
import six
from collections import OrderedDict
import sys, os
from math import sqrt, floor, ceil, log, exp, log2
from multiprocessing import Pool
import json

from fpylll import BKZ as BKZ_FPYLLL, GSO, IntegerMatrix, FPLLL
from fpylll.tools.quality import basis_quality
from fpylll.algorithms.bkz2 import BKZReduction

from custom_pump import my_pump_n_jump_bkz_tour

try:
	from g6k import Siever, SieverParams
	from g6k.algorithms.bkz import pump_n_jump_bkz_tour
	from g6k.utils.stats import dummy_tracer
except ImportError:
	raise ImportError("g6k not installed")


from keygen import NTRUKeyGenerator
from utils import parse_args, run_all, dump_basis, plot_gso, dump_blocksize, checkDSD, remove_zeros
import pre_processing
import post_processing

FPLLL.set_precision(120)


"""
Threshold on the length bound for the SRK event.
Interpreted as a*n^b
"""
thresholds = {
  "classic": [4,1],
  "classic_slice": [4,1],
  "phi": [4,1],
  "phi_projected": [4,1],
  "phi_projected": [4,3],
  "phi_projected_slice": [4,3],
  "phi_one": [4,1],
  "random": [4,1]
}


"""

# FIG.2:

for n in {101, 111, 121, 131, 141, 151, 161, 171, 181, 191, 201, 211}
python attack_ntru HRSS $n --bkz_alg=pump_n_jump --bkz_betas=64:110:1 --bkz_pre_beta=65 --sieve=bgj1 --threads=${nthreads} --verbose=True --dump=True --sieve_jump=2


#FIG. 3:
for n in {101, 111, 121, 131, 141, 151, 161, 171}
python attack_ntru HPS $n -q=512 --bkz_alg=pump_n_jump --bkz_betas=64:110:1 --bkz_pre_beta=65 --sieve=bgj1 --threads=${nthreads} --verbose=True --dump=True --sieve_jump=2

"""


class NTRUAttack:
	"""
	main NTRU attack class

	"""
	def __init__(self, params):
		self.n             = params['n']				# defines NTRU ring x^n - 1
		self.ntru_type     = params['ntru_type']		# HPS or HRSS
		self.q             = params['q']				# modulus (relevant for HPS, Challenge)
		self.lattice_type  = params['lattice_type']		# {classic, classic_slice,phi, phi_projected, phi_projected_slice, phi_one,challenge}
		self.nsamples      = params['nsamples']			# number of NTRU samples
		self.seed          = params['seed']				# NTRU key seed
		self.algbkz        = params['algbkz']			# {fpylll, pump_n_jump} = {enumeration, lattices}
		self.blocksizes    = params['blocksizes']		# range [a,b] with starting blocksize a, last blocksize b-1
		self.pre_blocksize = params['pre_blocksize']	# preproces lattice with enumeration up to blocksize pre_blocksize
		self.ntours        = params['tours']			# number of bkz tours

		self.nthreads      = params['threads']			# number of threads for sieving
		self.verbose       = params['verbose']			# verbose mode

		self.dump_basis    = params['dump']				# dump the basis after each blocksize together with the gso plot
		self.fileName      = params['filename']			# file to dump the basis and the gso (will be padded with the seed in the corresponding functions)

		#sieving related parameters

		self.sieve_alg       = params['sieve']			 # {nv, bgj1, gauss, hk3, bdgl}
		self.sieve_jump      = params['jump']			 # defines by how much we shift the context during a bkz tour
		self.pump_down_sieve = params['pump_down_sieve'] # flag whether to sieving in the pump down phase
		self.dim4free		 = params['dim4free']		 # [a,b] for the custom dim4free = a+b*blocksize

		self.pseudoinverse   = params['pseudoinverse']   # flag whether to use pseudoinverse of h (for hybrid)
		self.readfrom = params['readfrom']  			 # filename to read from a dumped basis

		self.c_in_p = params['check_in_preprocessing'] 	 # flag whether to check for SRK during the prepocessing or not

		if self.dim4free == None: self.dim4free=[11.5, 0.075] #defalut dim4free

		if params['h'] != None:
			self.h = json.loads(params['h'])  # for user-input h from NTRU Challenges
		else:
			self.h = None

		if self.pseudoinverse:
			self.lattice_type = "phi" #"phi_projected", "phi_projected_slice" might also work

		if self.readfrom == None:

			keyGenSuccess = False

			while not keyGenSuccess:
				self.generator  = NTRUKeyGenerator(self.ntru_type=="HRSS",self.n,self.pseudoinverse,self.q,self.seed,self.h)

				self.seedkeyGen = self.generator.newSeed()

				try:
					if self.lattice_type in ["classic", "classic_slice","challenge"]:
						self.basis = self.generator.getLattice(self.seedkeyGen)
					elif self.lattice_type in ["phi", "phi_projected", "phi_projected_slice"]:
						self.basis = self.generator.getLatticePhi(self.seedkeyGen)
					elif self.lattice_type == "random":
						self.basis = self.generator.getLatticeRandom()
					else: self.basis = self.generator.getLatticePhiOne(self.seedkeyGen)

					keyGenSuccess = True

				except ZeroDivisionError:
					self.seed += 512
					keyGenSuccess = False

			if self.lattice_type.endswith("_slice"):
				self.basis = pre_processing.sliceBasis(self.basis, params['sage'])

			if self.lattice_type=="challenge":
				self.basis = pre_processing.challenge(self.basis, params['sage'])

			if "projected" in self.lattice_type:
			  self.basis = pre_processing.projectAgainstOne(self.basis)

			self.B0 = IntegerMatrix.from_matrix(self.basis, int_type="long")

		else:

			self.B0 = IntegerMatrix.from_file("basis_dumps/"+self.readfrom)

		if self.B0.nrows <= 160:
			self.float_type = "long double"
		elif self.B0.nrows <= 450:
			self.float_type = "dd"
		else:
			self.float_type = "mpfr"

		self.M0 = GSO.Mat(self.B0, float_type=self.float_type)

		#LLL and remove lin. dependencies if present
		bkz = BKZReduction(self.M0)
		bkz.lll_obj() #
		B = remove_zeros(self.M0.B)


		M = GSO.Mat(B, float_type=self.float_type,
						U=IntegerMatrix.identity(B.nrows, int_type=B.int_type),
						UinvT=IntegerMatrix.identity(B.nrows, int_type=B.int_type))

		self.bkz = BKZReduction(M)

		if self.lattice_type=="challenge":
			self.sqnormH = sum([h_i*h_i for h_i in self.h])
			self.len_threshold = (2 * self.sqnormH * self.n) ** 2
		else:
			threshold_ = thresholds[self.lattice_type]
			self.len_threshold = (threshold_[0]*self.n**threshold_[1])

		self.param_sieve = SieverParams()
		self.param_sieve['threads'] = self.nthreads
		self.param_sieve['default_sieve'] = self.sieve_alg
		#self.param_sieve['saturation_ratio'] = 0.44
		self.g6k = Siever(M, self.param_sieve)
		print(self.g6k.params)

		if self.algbkz == "fpylll":
			self.M_ = self.bkz.M
		else:
			self.M_ = self.g6k.M

		self.B   = self.M_.B
		self.dim = self.B.nrows

		self.foundSubLattice = False		#
		self.blocksizeSubLattice = 0		# save blocksize when DSD happend
		self.blocksizeSecretKey = 0			# save blocksize when SKR happend


	def __call__(self):
		"""
		Calls lattice reduction.
		Prepocesses the basis if the basis was not dumped
		"""
		if self.readfrom == None:
			self.prebkz()
		if self.blocksizeSecretKey == 0:
			self.reduce()

	def prebkz(self):
		"""
		preprocess the lattice with enumeration
		"""
		if self.verbose:
			print('start preprocesssing the basis with prebeta = %d' % self.pre_blocksize)

		for b in range(7, self.pre_blocksize + 1):
			par = BKZ_FPYLLL.Param(
				b,
				strategies=BKZ_FPYLLL.DEFAULT_STRATEGY,
				max_loops=1,
				flags=BKZ_FPYLLL.MAX_LOOPS
			)
			self.bkz(par)

			if self.c_in_p:
				candidateShortest = self.getCandidateSKR()
				if self.checkSKR(candidateShortest):
					self.blocksizeSecretKey = b - 1

					candidateShortest = self.postProcessing(candidateShortest)

					if self.verbose:
						print('Found short vector at beta = %d' % (b - 1) )
						print(candidateShortest)
				break



	def reduce(self):
		"""
		Progressively reduce the basis.
		After each bkz tour check for DSD or SKR event
		"""

		T0 = time.time()

		for blocksize in self.blocksizes:

			if self.verbose: print("New round blocksize:", blocksize)

			if not self.foundSubLattice:
				foundSubLattice = checkDSD(self.M_)
				if foundSubLattice:  # exctract the dense sublattice and continue searching for a shortest vector in the extracted lattice
					self.B = self.B[:int(self.dim/2)]

					M_sublat = GSO.Mat(self.B, float_type=self.float_type,
									U=IntegerMatrix.identity(self.B.nrows, int_type=self.B.int_type),
									UinvT=IntegerMatrix.identity(self.B.nrows, int_type=self.B.int_type))

					self.bkz = BKZReduction(M_sublat)
					self.g6k = Siever(M_sublat, self.param_sieve)

					if self.algbkz == "fpylll":
						self.M_ = self.bkz.M
					else:
						self.M_ = self.g6k.M

					self.blocksizeSubLattice = blocksize - 1
					self.foundSubLattice = True

					if self.verbose:
						print('Found dense sublattice at beta = %d' % ( blocksize - 1) )

			candidateShortest = self.getCandidateSKR()

			if self.checkSKR(candidateShortest):
				self.blocksizeSecretKey = blocksize - 1

				candidateShortest = self.postProcessing(candidateShortest)

				if self.verbose:
					print('Found short vector at beta = %d' % (blocksize - 1) )
					print(candidateShortest)


				if not self.foundSubLattice:
					self.blocksizeSubLattice = -1

				if self.dump_basis:
					dump_blocksize( self.blocksizeSecretKey, self.blocksizeSubLattice, self.fileName, self.seed )

				break

			for t in range(self.ntours): # runs BKZ tours
				if self.algbkz == "fpylll":
					par = BKZ_FPYLLL.Param(blocksize,
										   strategies=BKZ_FPYLLL.DEFAULT_STRATEGY,
										   max_loops=8)
					self.bkz(par)

				elif self.algbkz == "pump_n_jump":
					my_pump_n_jump_bkz_tour(self.g6k, dummy_tracer, blocksize, jump=self.sieve_jump,
										 filename=self.fileName, seed=self.seed,
										 dim4free_fun="custom_dim4free_fun",
										 dim4free_param=self.dim4free,
										 extra_dim4free=0,
										 pump_params={'down_sieve': False},
										 goal_r0=self.len_threshold,
										 verbose=self.verbose)

			if self.dump_basis:
				dump_basis(self.M_.B, self.fileName+'_b_'+str(blocksize),seed=self.seed)
				plot_gso(self.M_, self.fileName+'_b_'+str(blocksize), seed=self.seed)

			if self.verbose:
				slope = basis_quality(self.M_)["/"]
				fmt = "{'lattype': '%10s',  'alg': '%15s', 'beta': %2d, 'slope': %.5f, 'total walltime': %.3f}" # noqa
				print(fmt % (self.lattice_type, self.algbkz + "+" + ("enum" if self.algbkz == "fpylll" else self.g6k.params.default_sieve),
								 blocksize, slope, time.time() - T0))


	def getCandidateSKR(self):
		"""
		Check whether lattices basis vectors are not among
		trivially short vectord (relevant for Coppersmith-Shamir lattice)
		"""
		d = int( self.M_.B.ncols/2 )
		badVectors = [
		2*d*[1],
		2*d*[-1],
		d*[1] + d*[0],
		d*[-1] + d*[0],
		d*[0] + d*[1],
		d*[0] + d*[-1]
		]
		i = 0
		while list(self.M_.B[i]) in badVectors or list( self.M_.B[i] )[:d] == [0]*d:
			i += 1

		return self.M_.B[i]

	def checkSKR(self, v):
		"""
		Check for the secret key recovery
		"""
		if self.verbose:
			print("Checking SKR", v)

		sqNorm = sum( v_i**2 for v_i in v )
		return sqNorm <= self.len_threshold
		#return v.norm() <= self.len_threshold


	def postProcessing(self, v):
		"""
		lattices of types {projected, challenge} require lifting
		"""
		if "projected" in self.lattice_type:
			v = post_processing.projected(v)
		if self.lattice_type=="challenge":
			v = post_processing.challenge(v, self.h, self.sqnormH)
		return list(v)


def main_call(params):
	attack_inst = NTRUAttack(params)
	return attack_inst()


if __name__ == '__main__':
	all_parsed_params = parse_args()
	run_all(main_call,all_parsed_params)
