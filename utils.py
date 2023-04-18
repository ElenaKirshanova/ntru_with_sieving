import argparse
import six
import copy
import sys, os
import re
from math import sqrt, floor, ceil, log, exp, log2
from datetime import datetime

from random import randint
from multiprocessing import Pool
import matplotlib.pyplot as plt

from fpylll.util import gaussian_heuristic
from fpylll import IntegerMatrix

import warnings



max_beta = 150 # maximal dimension we run SVP in
cross_over_beta = 70

def remove_zeros(B):
    """
    removes zero rows from matrix B
    """
    cr = 0
    for i in range(B.nrows):
      if not B[i].is_zero(): cr+=1

    B_ = [0]*cr
    cr = 0
    for i in range(B.nrows):
      if not B[i].is_zero():
        B_[cr] = B[i]
        cr+=1

    return IntegerMatrix.from_matrix(B_, int_type="long")

def testAndMakeDir(path):
  if not os.path.isdir(path):
      os.makedirs(path)

def checkDSD(M):
    """
        M: a GSO object
        checks if the DSD event happens
    """
    s_1 = 0
    s_2 = 0
    for i in range(M.B.nrows):
      norm = M.get_r(i,i)
      if i < M.B.nrows/2:
        s_1 += norm
      else:
        s_2 += norm

    return s_1 <= s_2

def dump_basis(B, filename, seed=None):
    """
    dumps integral basis B to a file named filename_seed#seedvalue.txt in the localpath/basis_dumps/
    """

    path = "basis_dumps/"
    testAndMakeDir(path)

    if not seed == None: filename += '_seed'+str(seed)+'.txt'
    else: filename += '.txt'

    d = B.nrows
    original_stdout = sys.stdout
    with open(path+filename, 'w') as f:
        sys.stdout = f
        #print('[')
        #for i in range(d):
        print(str(B))
        #print(']')

    sys.stdout = original_stdout
    return 1

def plot_gso(M, filename, seed=None):
    """
        dumps a plot of GSO norms of M to a file names filename_seed#seedvalue.txt in the
        localpath/gso_dumps/
    """

    d = M.B.nrows
    path = "gso_dumps/"
    testAndMakeDir(path)

    fig_name = filename
    if not seed == None: fig_name += '_seed'+str(seed)+'.png'
    else: fig_name += '.png'


    datapoints = []
    for i in range(d):
        val = M.get_r(i,i)

        if i < 3:
          print("GSO-DUMP", val)

        try:
            if val != 0:
                datapoints.append( log( M.get_r(i,i), 2 ) )
        except:
          print(val)

    fig, ax = plt.subplots( nrows=1, ncols=1 )
    ax.plot( datapoints )
    fig.savefig(path + fig_name)
    plt.close(fig)

    return 1

def dump_blocksize(skr, dsd, filename, seed):

  seed = seed - (seed % 10**8)

  path = "blocksize_dumps/"
  testAndMakeDir(path)

  filename += "_" + str(seed) + ".txt"

  with open(path+filename, "a+") as f:
    print( str(skr) + "\t" + str(dsd) + "\t" + datetime.now().strftime('%Y-%m-%d %H:%M:%S'), file = f )

def rough_estimate_on_betas(n,q):
    """
    use the output of this function in case the use did not provide us with blocksizes
    """
    if n<150: beta_low = 10
    else: beta_low = floor( 0.28*4.*n/( (log(q)/log(n))**2 + 1))
    return list(range(beta_low, max_beta))

def obtain_hrss_q_from_n(n):
    """
    HRSS modulus
    """
    return 2**ceil(7/2+log2(n))


def parse_args():
    parser = argparse.ArgumentParser(description='Parse NTRU attack params.')

    #main parameters
    parser.add_argument('ntru_type', type=str, help="HRSS or HPS or C")
    parser.add_argument('n', type=int,  help="ring dimension")
    parser.add_argument('-q', type=int, dest="q", default=None, help="NTRU modulus")
    parser.add_argument('--lattice_type', dest="lattice_type", default="phi_projected", help="{classic, classic_slice,phi, phi_projected, phi_projected_slice, phi_one, challenge}")
    parser.add_argument('--nsamples', type=int, default=None, dest="nsamples", help="Number of samples/rows of rot(h) used")
    parser.add_argument('--seed',  type=int, dest="seed", default=None, help="randomness seed")
    parser.add_argument('--pseudoinverse',   dest="pseudoinverse", default=False, help="If True, use pseudoinverse of h")
    parser.add_argument('--h', dest="h", default=None, help="Uses given input as h, instead of creating a random instance.")

    # number of runs, number of threads
    parser.add_argument('-t', '--trials', type=int, dest="trials", default=1,
                        help="number of experiments to run per dimension")
    parser.add_argument('-w', '--workers', type=int, dest="workers", default=1,
                        help="number of parallel experiments to run")
    parser.add_argument('--threads', type=int, dest="threads", default = 1, help="number of threads used by 1 worker")

    #bkz related params
    parser.add_argument('--bkz_alg', type=str, dest='algbkz', default="fpylll", help="{fpylll, pump_n_jump}")
    parser.add_argument('--bkz_betas', type=str, dest="blocksizes", default=None, help="bkz block sizes as string of the form: min_beta:max_beta:step")
    parser.add_argument('--bkz_pre_beta', type=int, dest="pre_blocksize", default=None, help="prepocessing block size")
    parser.add_argument('--bkz_tours', type=int, dest="tours", default=8, help="number of tours of bkz reduction")

    # sieving related params
    parser.add_argument('--sieve', type=str, dest="sieve", default="hk3", help="{nv, bgj1, gauss, hk3, bdgl}")
    parser.add_argument('--sieve_jump', type=int, dest="jump", default=1)
    parser.add_argument('--pump_down_sieve', dest = "pump_down_sieve", action='store_true')
    parser.add_argument('--dim4free', type=str, dest="dim4free", default=None, help="string of the form a:b. Dimensions for free are computed as a+b*blocksize (default: 11.5+0.075*blocksize)")
    # debug
    parser.add_argument('--verbose', dest="verbose", default=False, help="verbosity")
    parser.add_argument('--dry-run', dest="dry_run", default=False,
                        help="Show parameters that would be used but don't run any actual experiments.")
    parser.add_argument('--show-defaults', dest="show_defaults", action='store_true',
                        help="Show default parameters and exit.")

    # dump files
    parser.add_argument('--dump', dest='dump', default=False, help="flag to dump intermediate bases")
    parser.add_argument('--filename', dest='filename', default=None, help="prefix of the dump filenames")
    parser.add_argument('--readfrom', dest='readfrom', default=None, help="name of the file to read a basis from")

    # other
    parser.add_argument('--check_in_preprocessing', dest='check_in_preprocessing', default=False, help='check for SKR during the preprocessing phase or not')
    parser.add_argument('--sage', dest='sage', default='sage', help="command to run sage on your machine. only necessary, when using a sliced lattice")
    parser.add_argument('--hybrid', dest="hybrid", type=int, default=0, help='0 - no hybrid, 1 - plain hybrid, 2 - mitm hybrid. 1 and 2 require integer paramters r and k')
    parser.add_argument('--k', dest="k", type=int, default=0, help='number of guessed dimensions. Should be even for mitm hybrid')
    parser.add_argument('--r', dest="r", type=int, default=0, help='number of q-ary vectors we throw away in hybrid')

    args, unknown = parser.parse_known_args()


    fmt = "{key:%ds}: {value}"%20

    if len(unknown)>0:
        print('Parameters', unknown, 'are not recognized and will be ignored')

    all_defaults = {key: parser.get_default(key) for key in vars(args)}

    if args.show_defaults:
        for k, v in six.iteritems(all_defaults):
            print(fmt.format(key=k, value=v))
        exit(0)

    all_params = check_parsed_params(vars(args))

    if args.dry_run:
        for k, v in six.iteritems(all_params):
            print(fmt.format(key=k, value=v))
        exit(0)


    return all_params

def check_parsed_params(params):

    if not params['ntru_type'] in ["HRSS", "HPS", "C"]:
        raise ValueError("ntru_type=%s not recognized." % params['ntru_type'])

    if params['ntru_type']=="HPS" and params['q']==None:
        raise ValueError("Provide q for NTRU HPS!")
    elif params['ntru_type']=="HRSS": params['q'] = obtain_hrss_q_from_n(params['n'])

    if params['ntru_type']=="C" and params['q']==None:
        raise ValueError("Provide q for NTRU Challenge!")

    if params['ntru_type']=="C" and params['h']==None:
        raise ValueError("Provide h for NTRU Challenge!")

    if not params['lattice_type'] in ["classic", "classic_slice", "phi", "phi_projected", "phi_projected_slice", "phi_one", "random, challenge"]:
        raise ValueError("lattice_type=%s not recognized." % params['lattice_type'])

    if params['nsamples']==None: params['nsamples'] = params['n']
    else: assert(params['nsamples'] > 0 and params['nsamples']<=params['n'])

    if not params['algbkz'] in ["fpylll", "pump_n_jump"]:
        raise ValueError("algbkz=%s not recognized." % params['algbkz'])

    if params['pseudoinverse'] and params['ntru_type']=="HRSS":
        raise ValueError('pseudoinverse is implemented and tested only for HPS parameters')

    if params['blocksizes'] ==None:
        params['blocksizes'] = rough_estimate_on_betas(params['n'], params['q'])
    else: params['blocksizes'] = eval("range(%s)" % re.sub(":", ",", params['blocksizes']))

    if not params['dim4free'] ==  None:
        a, b = params['dim4free'].split(":")
        params['dim4free'] = float(a), float(b)

    assert(len(params['blocksizes'])>0)

    if params['pre_blocksize']==None:
        params['pre_blocksize'] = params['blocksizes'][0]-1

    if (params['dump'] and params['filename']==None):
        params['filename'] = 'n_'+str(params['n'])+'_lattype_'+str(params['lattice_type'])

    if (not params['readfrom']==None) and (not os.path.isfile("basis_dumps/"+params['readfrom'])):
         raise ValueError("file %s was not found." % params['readfrom'])

    if params['seed']==None:
        params['seed'] = randint(0, 2**64)

    # warn the caller is the choice of (bkz_alg, beta) is not optimal
    if (cross_over_beta in params['blocksizes']) and params['algbkz']=="fpylll":
        warnings.warn("Warning: enumeration in blocksizes" + str(params['blocksizes']) + 'will be slow. Run --bkz_alg="pump_n_jump" instead')

    if not int(params['hybrid']) in [0,1,2]:
        raise ValueError("unrecognized value for --hybrid. It can only be 0,1,2")
    if params['hybrid'] in [1,2] and params['k']==0:
        raise ValueError("provide --k parameter for the hybrid attack")
    if params['hybrid']==2 and not params['k']%2==0:
        raise ValueError("for mimt hybrid attack k should be even")


    return params

def run_all(f, params):
    jobs = []

    original_seed = params['seed']
    for t in range(params['trials']):
        params_  = copy.deepcopy(params)
        params_['seed'] = original_seed+t
        jobs.append(params_)
    if params['workers'] == 1:
        for job in jobs:
            res = f(copy.deepcopy(job))
    else:
        pool = Pool(params['workers'])
        pool.map(f, jobs)
