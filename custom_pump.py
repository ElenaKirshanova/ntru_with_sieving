"""

This is a modified pump.py file from
https://github.com/fplll/g6k/blob/master/g6k/algorithms/pump.py

The modifications include:
-- custom_dim4free_fun(): the user chooses their own [a,b] to have a+b*blocksize dim4free
-- we check for the DSD event inside the pump
-- in case of saturatin error we dump the GSO (for the purpose of understanding why this even happens frequently for ntru lattices)

"""
import sys, time
import six
from g6k.algorithms.workout import workout
from g6k.siever import SaturationError
import logging
from math import log

from utils import plot_gso, checkDSD

def print_pump_state(pump):
    pump.minl = min(pump.g6k.l, pump.minl)
    if pump.phase != "down":
        print("\r %3d: ↑%3d      " % (pump.r-pump.l, pump.g6k.r-pump.g6k.l), end=' ')
    else:
        print("\r %3d: ↑%3d ↓%3d " % (pump.r-pump.l, pump.r-pump.minl, pump.r-pump.g6k.l), end=' ')
    sys.stdout.flush()

def wrapped_sieve(pump, blocksize, filename, seed):
    if pump.phase == "init":
        alg = "gauss"
    else:
        alg = None

    cont = True
    try:
        with pump.g6k.temp_params(saturation_ratio=pump.g6k.params.saturation_ratio * pump.sat_factor):

            # Match lifting effort to insertion strategy
            pump.g6k(alg=alg, tracer=pump.tracer)

    except SaturationError as e:
        plot_gso(pump.g6k.M, filename+'_b_'+str(blocksize)+'_l'+str(pump.l)+'_r'+str(pump.r), seed=seed)
        if pump.saturation_error == "skip":
            pump.down_sieve = False
            logging.info("saturation issue: breaking pump.")
            cont = False
        elif pump.saturation_error == "weaken":
            logging.info("saturation issue: weakening pump.")
            pump.sat_factor /= 2.
        elif pump.saturation_error == "ignore":
            pass
        else:
            raise e

    if checkDSD(pump.g6k.M):
        print('DSD event inside pump')
        cont = False

    if pump.phase == "up" and (pump.max_up_time is not None):
        if pump.max_up_time < time.time() - pump.up_time_start:
            cont = False

    if pump.goal_r0 is not None:
        pump.g6k.insert_best_lift(scoring_goal_r0, aux=pump)

        if (pump.g6k.M.get_r(pump.kappa, pump.kappa) <= pump.goal_r0):
            cont = False

    return cont


def dim4free_wrapper(dim4free_fun, dim4free_param, blocksize):
    """
    Deals with correct dim4free choices for edge cases when non default
    function is chosen.

    :param dim4free_fun: the function for choosing the amount of dim4free
    :param blocksize: the BKZ blocksize

    """
    if blocksize < 40:
        return 0
    dim4free = dim4free_fun(blocksize, dim4free_param)
    return int(min((blocksize - 40)/2, dim4free))


def default_dim4free_fun(blocksize):
    """
    Return expected number of dimensions for free, from exact-SVP experiments.

    :param blocksize: the BKZ blocksize

    """
    return int(11.5 + 0.075*blocksize)

def custom_dim4free_fun(blocksize, dim4free_param):
    """
    Return expected number of dimensions for free, from exact-SVP experiments.

    :param blocksize: the BKZ blocksize

    """
    return int(dim4free_param[0] + dim4free_param[1]*blocksize)

def scoring_(i, nlen, olen, aux):
    return i == aux.kappa and nlen < aux.goal_r0


def scoring_down(i, nlen, olen, aux):
    if i < aux.insert_left_bound or nlen >= olen:
        return False
    return log(olen / nlen) - i * log(aux.prefer_left_insert)

def pump(g6k, tracer, kappa, blocksize, dim4free, filename, seed, down_sieve=False,                 # Main parameters
         goal_r0=None, max_up_time=None, down_stop=None, start_up_n=30, saturation_error="weaken",  # Flow control of the pump
         increasing_insert_index=True, prefer_left_insert=1.04,                                     # Insertion policy
         verbose=False,                                                                             # Misc
         ):
    """
    Run the pump algorithm.

    :param g6k: The g6k object to work with
    :param tracer: A tracer for g6k
    :param kappa: beginning of the block
    :param blocksize: dimension of the block (r=kappa+blocksize)
    :param dim4free: number of ``dimension for free'' [Ducas, Eurcrypt 2018]: Sieve context [l,r] where l=kappa+dim4free
    :param down_sieve: re-sieve after each insert during the pump-down phase.  (stronger reduction,
        slower running time)
    :param goal_r0: an extra hook to always insert at position kappa if this goal length can be met
        by a lift.  Quit when this is reached.
    :param max_up_time: For balancing BKZ time with SVP call time in LWE.  Stop pumping up when this
        time has elapsed.
    :param down_stop: stop inserts during pumping down after index kappa+down_stop (to control
        overheads of insert in large dimensional lattices)
    :param start_up_n: Initial sieve-context dimension for pumping up (starting at 1 incurs useless overheads)
    :param saturation_error: determines the behavior of pump when encountering saturation issue {"weaken",
        "skip", "ignore", "forward"}
    :param increasing_insert_index: During pump-down, always insert on the right side of the previous insertion.
    :param prefer_left_insert: Parameter theta from the paper (Sec 4.4) for scoring insertion candidates.
    :param verbose: print pump steps on the standard output.

    """
    pump.r = kappa+blocksize
    pump.l = kappa+dim4free  # noqa

    g6k.shrink_db(0)
    g6k.lll(kappa, pump.r)
    g6k.initialize_local(kappa, max(pump.r-start_up_n, pump.l+1), pump.r)

    pump.sat_factor = 1.
    pump.up_time_start = time.time()
    pump.insert_left_bound = kappa
    pump.minl = g6k.l

    for key in ('kappa', 'down_sieve', 'goal_r0', 'g6k', 'tracer',
                'max_up_time', 'saturation_error', 'verbose', 'prefer_left_insert'):
        setattr(pump, key, locals()[key])

    if down_stop is None:
        down_stop = dim4free

    with tracer.context(("pump", "beta:%d f:%d" % (blocksize, dim4free))):
        with g6k.temp_params(reserved_n=pump.r-pump.l):
            pump.phase = "init"
            wrapped_sieve(pump, blocksize,filename, seed)  # The first initializing Sieve should always be Gauss to avoid rank-loss

            pump.phase = "up"
            # Pump Up
            while (g6k.l > pump.l):
                g6k.extend_left(1)

                if verbose:
                    print_pump_state(pump)
                if not wrapped_sieve(pump, blocksize,filename, seed):
                    break

            if goal_r0 is not None and (g6k.M.get_r(kappa, kappa) <= goal_r0):
                return

            # Pump Down
            pump.phase = "down"
            while (g6k.n > 1) and (pump.insert_left_bound <= kappa+down_stop):
                # (try to) Insert
                ii = g6k.insert_best_lift(scoring_down, aux=pump)
                if ii is not None and increasing_insert_index:
                    pump.insert_left_bound = ii + 1

                else:
                    g6k.shrink_left(1)

                if goal_r0 is not None and (g6k.M.get_r(kappa, kappa) <= goal_r0):
                    break

                # Sieve (or Shrink db)
                if verbose:
                    print_pump_state(pump)
                if not pump.down_sieve:
                    g6k.resize_db(max(500, g6k.db_size() / g6k.params.db_size_base))
                elif not wrapped_sieve(pump, blocksize,filename, seed):
                    break

def my_pump_n_jump_bkz_tour(g6k, tracer, blocksize, jump=1,
                         filename=None, seed=None,
                         dim4free_fun=default_dim4free_fun, dim4free_param=[11.5, 0.075], extra_dim4free=0,
                         pump_params=None, goal_r0=0., verbose=False):
    """
    Run a PumpNjump BKZ-tour: call Pump consecutively on every (jth) block.

    :param g6k: The g6k object to work with
    :param tracer: A tracer for g6k
    :param blocksize: dimension of the blocks
    :param jump: only call the pump every j blocks
    :param dim4free_fun: number of dimension for free as a function of beta (function, or string
        e.g. `lambda x: 11.5+0.075*x`)
    :param extra_dim4free: increase the number of dims 4 free (blocksize is increased, but not sieve
        dimension)
    :param pump_params: parameters to pass to the pump
    """
    if pump_params is None:
        pump_params = {"down_sieve": False}

    if "dim4free" in pump_params:
        raise ValueError("In pump_n_jump_bkz, you should choose dim4free via dim4free_fun.")

    d = g6k.full_n
    g6k.shrink_db(0)
    g6k.lll(0, d)
    g6k.update_gso(0, d)


    if isinstance(dim4free_fun, six.string_types):
        dim4free_fun = eval(dim4free_fun)

    dim4free = dim4free_wrapper(dim4free_fun, dim4free_param, blocksize) + extra_dim4free
    blocksize += extra_dim4free

    indices  = [(0, blocksize - dim4free + i, i) for i in range(0, dim4free, jump)]
    indices += [(i, blocksize, dim4free) for i in range(0, d - blocksize, jump)]
    indices += [(d - blocksize + i, blocksize - i, dim4free - i) for i in range(0, dim4free, jump)]

    pump_params["down_stop"] = dim4free + 3

    for (kappa, beta, f) in indices:
        if verbose:
            print("\r k:%d, b:%d, f:%d " % (kappa, beta, f), end=' ')
            sys.stdout.flush()

        pump(g6k, tracer, kappa, beta, f, filename, seed, **pump_params)
        g6k.lll(0, d)
        if g6k.M.get_r(0, 0) <= goal_r0:
            return

    if verbose:
        print("\r k:%d, b:%d, f:%d " % (d-(blocksize-dim4free), blocksize-dim4free, 0), end=' ')
        sys.stdout.flush()
    pump_params["down_stop"] = blocksize-dim4free
    pump(g6k, tracer, d-(blocksize-dim4free), blocksize-dim4free, 0,  filename, seed, **pump_params)
    if verbose:
        print()
