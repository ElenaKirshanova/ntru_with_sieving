# ntru_solver
Solving NTRU with lattice reduction.
This repository contains the scripts accompanying the article

**Attacking NTRU in Practice: Cyclotomic Ring, Almost-Parallel Hints, and Sieving**


# Requirements

* [fpylll](https://github.com/fplll/fpylll)
* [g6k](https://github.com/fplll/g6k)
* [SageMath 9.3+](https://www.sagemath.org/) (required only if one wants to use either of the lattices types {classic_slice,phi_projected_slice, challenge} or to solve [NTRU Challenges](https://web.archive.org/web/20160310141551/https://www.securityinnovation.com/uploads/ntru-challenge-parameter-sets-and-public-keys-new.pdf))


# Description of files
Short description of the content:
* `attack_ntru.py` main file to run lattice reduction attack on NTRU
* `keygen.py` generates NTRU-HRSS and NTRU-HPS parameters
* `custom_pump.py` is a modified pump function from the  [g6k](https://github.com/fplll/g6k) library
* `pre_processing.py` preprocesses lattice bases of type 'slice', 'projected', 'challenge'
* `post_processing.py` post-processes the shortest vector found during the reduction
* `utils.py` contains helping functions
* `dual_basis.sage` computes the dual basis. Scales the basis such that it is integral
* folder `gso_dumps` contains GSOs of dumped bases during the reduction for HPS and HRSS experiments (see Figures 1 and 2 in the accompanying paper). The file names are of the form
`n_${n}_${lattype}_b_${beta}_seed${seed}`, e.g. `n_121_lattype_phi_projected_b_2_seed8494989096862686174`.
* `seeds.txt` contains seeds to recreate Figures 1 and 2
* `fit.sage ` interpolates the lines for Figures 1 and 2
* `challenge181.sage` verifies the solution for the [NTRU-181 Challenges](https://web.archive.org/web/20160310141551/https://www.securityinnovation.com/uploads/ntru-challenge-parameter-sets-and-public-keys-new.pdf)) that we found


# How to use

Run `attack_ntru.py` with the following parameters

* {HRSS, HPS, C} : NTRU type as per [NIST 3rd Round Specification document](https://ntru.org/f/ntru-20190330.pdf) or as per NTRU Challenge from [Securty Innovation Inc.](https://web.archive.org/web/20160310141551/https://www.securityinnovation.com/uploads/ntru-challenge-parameter-sets-and-public-keys-new.pdf) Obligatory parameter
* `n` : defines NTRU ring x^n - 1. Integer. Obligatory parameter
* `-q`: NTRU modulus
* `--lattice_type` :  either of the following lattices {classic, classic_slice,phi, phi_projected, phi_projected_slice, challenge}
* `--nsamples` : Number of samples/rows of rot(h) used. Integer
* `--seed`: randomness seed to the NTRU key generation.
* `----pseudoinverse`: If True, use pseudoinverse of h as described [here](https://csrc.nist.gov/CSRC/media/Events/third-pqc-standardization-conference/documents/accepted-papers/nguyen-boosting-hybridboost-pqc2021.pdf)

* `-t` (or `--trials`): number of experiments to run per NTRU dimension (same NTRU type)
* `-w` (or `--workders`): number of parallel experiments to run
* `--threads`: number of threads used by 1 worker

BKZ-related parameters:

* `--bkz_alg`: either of following {fpylll, pump_n_jump}. Choose fpylll for enumeration, pump_n_jump for sieving as SVP oracles
* `--bkz_betas`: string of the form bkz_beta_start:bkz_beta_end:beta_step
* `--bkz_pre_beta`: blocksize for preprocessing (run BKZ with enumeration up to bkz_pre_beta (inclusive)). Integer
* `--bkz_tours`: number of bkz tours per blocksize. Integer


Sieving related parameters:

* `--sieve`: either of the following sieving algorithms {nv, bgj1, gauss, hk3, bdgl}
* `--sieve_jump`: determines by how much dimension we move between the sieving context. Integer.
* `--pump_down_sieve`: boolean that determines whether we sieving during the pump down
* `--dim4free`: string of the form a:b that determines the dimensions for free as a+b*blocksize.

* `--dump`: boolean that determines whether to dump the current basis after each blocksize is finished. Dumps bases to folder basis_dumps/
* `--filename`: name of the file to dump bases. The file will be named $filename+seed.txt
* `--readfrom`: filename to read the input basis from. It will be read from basis_dumps/$readfrom

* `--check_in_preprocessing`: boolean that determines whether to check for SKR during the preprocessing
* `--sage`: command to run sage on your machine. only necessary, when using a sliced lattice


# Experiments

To attack an NTRU-HRSS instance over x^111 - 1 using enumeration simply run
```
python attack_ntru.py HRSS 111 --verbose=True
```

It takes approximately half a minute on a laptop.
The script generates a random instance of NTRU-HRSS parameters according to the [NIST 3rd Round Specification document](https://ntru.org/f/ntru-20190330.pdf).

To generate an NTRU instance with a specific seed run

```
python attack_ntru.py HRSS 111 --seed=1
```

This run takes 34sec. on Apple M1 Pro and the shortest vector is found at blocksize 11.


To attack an NTRU-HPS instance over x^121 - 1 with q=512 run
```
python attack_ntru.py HPS 121 -q=512 --alg_bkz=pump_n_jump --verbose=True
```
It takes approximately 1h30 on a laptop to find a shortest vector.

To recreate the black dots from Figure 1 (NTRU-HRSS), we run the following command (for each n from the set {101, 111, 121, 131, 141, 151, 161, 171}) with 20 parallel runs starting with seed=14498356423527637276, e.g. for n=151:
```
python attack_ntru.py HRSS 151 --lattice_type='classic' --bkz_betas=20:50:1 --seed=14498356423527637276 --verbose=True --dump=True --workers=20 --trials=32
```

For larger n's, e.g. n=171, we used sieving (hense --alg_bkz=pump_n_jump) in parallel (--threads=10) for betas>=65 and for preprocessing with enumeration for smaller beta:
```
python attack_ntru.py HRSS 171 --lattice_type='classic' --bkz_pre_beta=64 --bkz_betas=65:70:1 --seed=14498356423527637276 --verbose=True --dump=True --threads=10 --workers=20 --trials=32
```

To recreate blue squares (phi_projected lattice) we run for e.g. n=191:
```
python attack_ntru.py HRSS 191 --lattice_type='phi_projected' --bkz_pre_beta=64 --bkz_betas=65:90:1 --seed=2584232668076212209 --verbose=True --dump=True --threads=20 --workers=10 --trials=20
```

To experiments on HPS parameters one needs to provide the modulus q as well, e.g.
```
python attack_ntru.py HPS 161 -q=512 --lattice_type='phi_projected' --bkz_pre_beta=64 --bkz_betas=65:110:1 --seed=8550301121086457479 --verbose=True --dump=True --threads=20 --workers=10 --trials=32
```
recreates the blue square corresponding to n=161 in Figure 3.

The seeds for all dimensions can be found in seeds.txt

To attack NTRU Challenges from [Security Innovation Inc.](https://web.archive.org/web/20160310141551/https://www.securityinnovation.com/uploads/ntru-challenge-parameter-sets-and-public-keys-new.pdf) add NTRU type `C`, `--h` option to include the public key provided by the challenges, and `--lattice_type=challenge` like so
```
python attack_ntru.py C 113 -q=1024 --verbose=True --threads=256 --bkz_alg=pump_n_jump  --lattice_type=challenge --dump=True --h="[1021, 344, 107, 401, 913, 817, 980, 271, 512, 110, 957, 736, 920, 275, 807, 210, 41, 65, 280, 161, 588, 561, 312, 252, 472, 158, 999, 963, 906, 519, 40, 447, 885, 970, 1015, 234, 534, 512, 47, 899, 630, 654, 548, 113, 410, 798, 908, 758, 443, 44, 690, 380, 672, 143, 732, 142, 191, 129, 97, 606, 993, 944, 399, 138, 260, 745, 561, 428, 813, 459, 853, 604, 990, 813, 533, 37, 334, 143, 601, 389, 473, 713, 770, 7, 641, 246, 417, 782, 731, 562, 908, 789, 285, 320, 534, 271, 676, 1, 813, 788, 462, 169, 634, 986, 592, 168, 110, 648, 310, 502, 946, 1001, 772]"
```

# Solving the 181-NTRU Challenge


The solution to the 181 Challenge if the form (g, f) is
```
[-1, 0, -1, -1, -1, -1, 0, 1, -1, 1, 1, 0, -1, -1, 0, 1, 0, 1, 1, -1, 1, 1, 0, 1, -1, 1, 1, 0, 0, 0, -1, 1, 0, -1, 0, -1, -1, -1, 0, -1, 0, -1, 0, 1, 1, 1, 1, -1, -1, 1, -1, 0, -1, -1, 0, 1, 1, 1, 1, -1, -1, 0, 0, 0, -1, -1, 0, 0, -1, 1,
1, 0, 1, 0, -1, 0, 1, 0, 0, 1, 1, -1, -1, 1, -1, 1, 0, 0, 0, -1, 0, 0, -1, -1, 1, 1, -1, 0, 0, 0, 1, -1, 1, 0, -1, 1, 1, -1, -1, 0, -1, -1, 0, -1, 0, 1, -1, 1, 0, 0, 1, -1, 0, 1, 0, -1, -1, 0, 1, 0, 1, -1, 1, -1, -1, -1, 0, 1, 1, -1, 1, 0, 0, -1, -1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, -1, 1, 1, 1, -1, -1, 0, -1, 1, 0, 0, 0, 1, 1, 1, -1, -1, 1, -1, 0, -1, 1, 1, 1, 0, 1, 1, 0, -1, -1, 0, 0, -1, 0, 0, -1, -1, 0, -1, -1, -1, 1, -1, 0, -1, -1, -2, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 2, 0, 0, 1, 2, 1, 0, -1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 2, 2, 1, 0, 0, 0, 1, 2, 3, 0, 0, 1, 0, 0, 0, 0, 1, -2, -2, 0, -1, 0, 0, 0, 0, -2, 1, 1, -1, 1, 0, 0, 0, -1, -1, 0, -1, -2, -1, 0, 0, 0,
0, 0, -1, -1, 0, 1, 0, 0, 0, 1, -1, -1, -1, 2, -2, -1, 0, 0, -1, -1, -1, 1, -2, 1, -1, 0, 0, 0, 1, 1, -1, 0, 0, 0, 1, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 0, 1, 1, -1, 0, 0, 2, 1, 1, 2, -1, 1, 0, 0, 1, 0, -1, 1, 0, 0, 0, 0, 0, 0, -1]
```
