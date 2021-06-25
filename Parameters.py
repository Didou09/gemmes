# -*- coding: utf-8 -*-
"""
PARAMETERS FOR GOODWIN-KEEN LIKE MODELS
"""


import numpy as np
import FunctionsGoodwin as FG


_PARAMSET = 'v0'
_FLATTEN = True


# #############################################################################
# #############################################################################
#                       Default dict of parameters 
# #############################################################################


# Default parameter dict, inc. all parameters
_DEF_PARAM = {

    'economy': {

        # SET OF EQUATIONS ########################################################
        'EquationSet': FG.f_reduced,    # The system of equation to solve

        ### Numerical PARAMETERS ##################################################
        # Time Discretisation ###
        'Tmax': 100,                    # Time simulated (in years)
        'time_safety': 30,              # Factor of time increment reduction. The system calculate an increment dtmax, then divide it by this value

        # Space discretisation/Number of systems 
        'Nx': 2,                        # Number of parrallel economy (that can be coupled)

        # DATA STORAGE #####
        'StorageMode': 'Full',          # Keep all iterations of only few
        'Tstore': 0.1,                  # Time between two storage (if StorageMode=full, it goes to dt)
        'Save': False,
        'printblabla': True,

        ### EACH MODULES AND VALUES ###############################################
        # Population evolution 
        'beta': 0.025,                  # Rate of population growth     (y^-1)
        'PopSat': 10000,                # Ratio between initial and maximal population 
        'alpha': 0.02,                  # Rate of productivity increase (y^-1)

        # Capital properties
        'delta2': 0.005,                # Rate of ALL capital depletion     (y^-1)
        'delta1': 0.005,                # Rate of WORKING capital depletion (y^-1)
        'pK': 0,                        # Price consumption of capital

        'nu': 3,                        # Kapital to output ratio as in Leontiev. CAREFUL IN CES this is 1/A
        'eta': 10000,                   # CES Only 1/(1+substituability)
        'b': .5,                        # CES parameter : capital part of the production 
        'z': 1,                         # Markup on salary estimation

        # INTEREST / Price
        'r': .03,                       # Interest at the bank
        'eta': .192,                    # Typical rate for inflation
        'mu': 1.3,                      # Mark-up of price
        'gamma': 1,                     # Money-illusion

        ### FUNCTIONS AND THEIR PARAMETERS ########################################
        # PHILIPS CURVE (employement-salary increase)
        'phiType': 'divergent',         # 'divergent' , 'linear' modes
        'phinul': 0.04,                 # Unemployement rate at which there is no salary increase with no inflation
        'phislope': 2,                  # slope in a linear model (power in a diverging model)         

        # KEEN INVESTMENT FUNCTION (profit-investment function)
        'k0': -0.0065,
        'k1': np.exp(-5),
        'k2': 20,

        # LINEAR DIVIDENT PROFITS 
        'div0': 0.138,                  # Part of GDP as dividends when pi=0
        'div1': 0.473,                  # Slope
        },

    'climate': {

        # Coupling and climate ####################################################
        # Coupling Effets
        'g1': .0,                       # GLOBAL         EFFECTS OF LAMBDA
        'g2': .00,                      # WITH NEIGHBORS EFFECTS OF LAMBDA
        'muI': 0,                       #
        'muN': 0,                       #

        # Recrutment-Firing dynamics
        #tau: np.linspace(0.8,1.1,p['Nx'])
        'tauR': 0.0, #85                # Typical time for recruitement (y)
        'tauF': 0.0,                    # Typical time for firing       (y)

        'tauL' : 0.0,                   # Typical relaxation time on employement perception
        'tauLam': 2.,
        'tauK' : 2,                     # Typical relaxation time on

        ### INITIAL CONIDITONS
        'CO2at_ini': 851,
        'CO2up_ini': 460,
        'CO2lo_ini': 1740,
        'T_ini': 0,

        # Climate model
        'Phi12' : .024,                 # Transfer of carbon from atmosphere to biosphere
        'Phi23' : .001,                 # Transfer from biosphere to stock

        'C'   : 1/.098,                 # Heat capacity of fast-paced climate
        'C0'  : 3.52,                   # Heat capacity of inertial component of climate
        'gamma': 0.0176,                # Heat exchange coef between layer !!CONFLICT?!
        'Tsens': 3.1,                   # Climate sensibility

        'dFexo': 0,
        'FexoMax': 0,
        'F2CO2' : 3.681,                # W/M2, doubling CO2 impact on forced radiations

        # CLIMATE-INTEGRATION
        '''
        Parameters from Debt&Damages
        '''

        'theta': 2.6,                   # Convexity on abattement cost function
        'dsigma': - 0.001,              # Variation rate of the growth of emission intensity
        'dPBS' : - 0.005,               # Growth rate of back-stop technology price
        'dPc' : 0, #[CORRECT]           # Growth rate of back-stop technology price
        'dEland': - 0.022,              # Growth rate of land use change in CO2 emission

        # Damage function (on GDP)
        '''
        D: 1 - (1 + param['pi1']*T + param['pi2']*T**2 + param['pi3']*T**param['zeta'] )**(-1)
        '''
        'pi1': 0,                       # Linear temperature impact
        'pi2': .00236,                  # Quadratic temperature impact
        'pi3': .00000507,               # Weitzmann Damage temperature impact
        'zeta': 6.754,                 # Weitzmann impact
        'fk': 1/3,                      # Fraction of environmental damage allocated to the stock of capital
    },
}


# ####################################################
# Creat a dict of different versions of parameter sets

_DPARAM = {
    'v0': {k0: dict(v0) for k0, v0 in _DEF_PARAM.items() if k0 in ['economy']},
    'v1': {k0: dict(v0) for k0, v0 in _DEF_PARAM.items()},
    'GreatAuthor2019': {
        k0: dict(v0) for k0, v0 in _DEF_PARAM.items() if k0 in ['economy']
    },
}


_DPARAM['GreatAuthor2019']['economy']['b'] = 0.
_DPARAM['GreatAuthor2019']['economy']['eta'] = 0.192


"""
f_reduced         #WORKING
f_ces_reduced

f_full_leontiev
f_full_CES

f_reduced_relax

f_GEMMES
f_GEMMES_CLIM
f_GEMMES_CORE
"""

# #############################################################################
# #############################################################################
#                       Utilities 
# #############################################################################


def _check_inputs(paramset=None, flatten=None):

    # paramset
    if paramset is None:
        paramset = _PARAMSET
    c0 = isinstance(paramset, str) and paramset in _DPARAM.keys()
    if not c0:
        ls = ['\t- {}'.format(kk) for kk in sorted(_DPARAM.keys())]
        msg = (
            "Arg paramset must be a valid predefined parameter set!\n"
            + "\n".join(ls)
            + "\nYou provided: {}".format(paramset)
        )
        raise Exception(msg)

    # flatten
    if flatten is None:
        flatten = _FLATTEN
    c0 = isinstance(flatten, bool)
    if not c0:
        msg = (
            "Arg flatten must be a bool!\n"
            + "\nYou provided: {}".format(flatten)
        )
        raise Exception(msg)

    return paramset, flatten


# #############################################################################
# #############################################################################
#           Choose which version of the dict of parameters to use 
# #############################################################################


def initparams(paramset=None, flatten=None):
    """
    Create a dictionnary containing all the parameters necessary for simulation
    Their description is in comments.

    parameters
    ----------
    paramset:   None / str
        Flag indicating which predefined set of parameter to pick
        Defaults to 'v0'
    flatten:    None / bool
        Flag indicating whether to flatten the param dict
        Used for retro-compatibility
        Default to True

    """

    # ------------
    # Check inputs
    paramset, flatten = _check_inputs(
        paramset=paramset, flatten=flatten,
    )

    # ------------
    # Dictionnary of parameters (copy to avoid modifying the original)
    param = {k0: dict(v0) for k0, v0 in _DPARAM[paramset].items()}

    # ------------
    # return in desired format
    if flatten is True:
        # flatten dict of parameters for retro-compatibility
        lk0 = list(param.keys())
        paramflat = dict(param[lk0[0]])
        for k0 in lk0[1:]:
            paramflat.update(param[k0])
        return paramflat
    else:
        return param
