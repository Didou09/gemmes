#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
The code is importing some other libraries so that we can use their functions.
"""

# Standard library
import os
import sys
import argparse


# common
import time                     # Time (run speed) printing
import numpy as np              # Array/Numerical operations
import matplotlib.pyplot as plt # Plot elements 
# plt.close('all')                # Remove to 


# gemmes-specific
import Parameters as Par      # All the parameters of the system
import FunctionsGoodwin as FG # Core of models
import plots as plts          # Already written plot functions
import Miscfunc as M          # All miscellaneous functions


_ROOTFOLD = os.path.dirname(__file__)
_PLOT = True


### WELCOME MESSAGE ###########################################################
"""
GOODWIN-TYPE RESOLUTION ALGORITHM

This python code has been created to simulate Goodwin-Keen simulations, with all the variants to it.
Do not worry about the number of section : the architecture is thought to be user-friendly and coder-friendly

WHAT ARE THE SPECIFICITIES OF THIS CODE?
* The algorithm solve p['Nx'] system in parrallel : the idea is that later we will couple the systems altogether ("spatial/network properties").
* For now, there is no much coupling, but you can use it to test many initial conditions and array of parameters values
* The steps are :
    - Creation of the parameters dictionnary p ( in Parameters.py for initial values)
    - Creation of initial conditions in the dictionnary ic
    - Translation into machine friendly variable y
    - Calculation of all the dynamics variable Y_s (the temporal loop)
    - Translation into user-friendly result in the dictionnary r
    - Additional analysis (period, slow enveloppe mouvement...)
    - Plot of variables and analysis through r

HOW TO IMPLEMENT YOUR TOY-MODEL?
 * Add the parameters in Parameters.py (check if it isn't already, or if the name isn't already taken)
 * Add the temporal component in FG.preparedT for dt selection
 * Add the related quantity in the dictionnary of initial condition ic
 * Add a section to prepareY
 * Add the equations in a f_THENAMEOFTHEFUNCTION
 * Add the redevelopment in FG.expand
 * have fun !

WHAT I AM (Paul) LOOKING FOR IN FURTHER DEVELOPMENT
* Check the functions in FunctionGoodwin, I've done some mistakes in more complex models (typically coping with collapse)
* Add spatial operators which are non unstable ( I have implicit scheme elsewhere but it's a different kind of resolution, and much more work when you change a model)
* Have a list of all models existing in this framework (copingwithcollapse, Harmoney, predatory-prey...) and code them
* The plots are UGLY ! Let's do something nicer
* As Iloveclim and Dymends are not in python, prepare some bindings
"""


# #############################################################################
# #############################################################################
#                   Ulities
# #############################################################################


def _str2bool(arg):
    if isinstance(arg, bool):
        return arg
    elif arg.lower() in ['yes', 'true', 'y', 't', '1']:
        return True
    elif arg.lower() in ['no', 'false', 'n', 'f', '0']:
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected!')


def _check_inputs(
    plot=None,
    rootfold=None,
):

    if plot is None:
        plot = _PLOT
    if rootfold is None:
        rootfold = _ROOTFOLD

    return plot, rootfold



# #############################################################################
# #############################################################################
#                   Main function
# #############################################################################


def run(
    plot=None,
    rootfold=None,
):
    """ Run the simulation

    Details

    parameters
    ----------
    plot:       bool
        Flag indicating whether to plot figure
    rootfold:   str
        Path where to save data


    """

    # -----------------
    #  check inputs

    plot, rootfold = _check_inputs(
        plot=plot,
        rootfold=rootfold,
    )

    # -----------------
    #  go on

    print('List of existing set of equations  :')
    for f in [ f for f in dir(FG) if 'f_' in f] : print(f)
    print(3*'\n')


    # -----------------
    # PARAMETERS INITIALISATION
    """
    This part create 'param' (parameters dic) and 'op' (operators dic) depending of the system size and properties
    Typically, param contains the "default parameters" values and should be open in another windows,
    I edit them after param initialisation if I want to do some original runs
    """

    param = Par.initparams()              # PARAMETERS IN THIS FILE ( PARAMETERS )

    param['lambdamax']  = 1-10**(-2)          # Variation of lambda on the simulation 

    # The good place to add a loop on certain parameters if you need 
    param = M.preparePhilips(param)    # Determining Philips curve and Solow point
    param = M.preparedT(param)         # Determine the best value for dt and associates
    op= []#M.prepareOperators(p)  # Spatial operators initialisation

    # -----------------
    # initial conditions
    """
    The initial state vector at t=0 of the system in a "readable" language.
    Depending of the equation set you use, all of them will not necessary be used.
    (example, a reduced set will read ic['d'], a complete set will read ic['D'])
    """

    v1 = np.ones(param['Nx']) ### Pratical notation for more readable code : an array with 1 everywhere

    # Dictionnary of initial conditions
    ic = {
        'a': 1.*v1,                             # Initial productivity per worker
        'N': 1.*v1,                             # Initial quantity of people
        'd': 0*1.53*v1,
        'omega': 1.*v1*param['omega0'],         #*np.arange(0,1,p['Nx'])
        'lambda': np.linspace(                  # Distribution of initial lambda
            param['lambdamin'], param['lambdamax'], param['Nx'],
        ),
    }

    ic['L'] = ic['omega']*ic['N']
    ic['K'] = ic['a']*ic['L']*param['nu']               # Beginning at equilibrium
    ic['W'] = param['omega0'] *ic['K']/param['nu']      # Initial distribution of salary

    '''
    ic['D']     =1.*v1            # Debt
    ic['sigma'] =1.*v1            # Intensity of emission in production
    ic['P']     =1.*v1            # Global price
    ic['Pbs']   =547.22*v1            # Price Backstop
    ic['Pc']    =1.*v1            # Carbon price


    ic['Fexo']  =1.*v1            # Exogenous forcing intensity
    ic['Eland'] =2.6*v1           # Emission from land
    ic['CO2at'] =1.*v1*param['CO2at_ini']           # CO2 in atmosphere
    ic['CO2up'] =1.*v1*param['CO2up_ini']           # CO2 in upper ocean and biosphere
    ic['CO2lo'] =1.*v1*param['CO2lo_ini']          # CO2 in deep ocean
    ic['T']     =1.*v1           # Atmospheric temperature
    ic['T0']    =1.*v1           # Deeper ocean temperature
    '''
    ic['t']     = 1.*v1           # Year in the simulation 

    y, param = FG.prepareY(ic, param) ### The vector y containing all the state of a time t.

    # -----------------
    # Initialisation of data storage
    """
    This part is very close to the numerical resolution core. No need to delve in it
    It looks "complex" because it allows the record of few timestep only
    """
    Y_s = np.zeros((param['Nvar'], param['Ns'], param['Nx']))     # stock dynamics
    Y_s[:, 0, :] = 1*y                                 # first stock
    t            = 0
    t_s          = np.zeros(param['Ns'])                   # stock time
    tprevious    = 0                                   # deltatime between two records
    idx          = 0                                   # index of the current record

    # -----------------
    # Simulation iteration
    """
    The RK4 core, with partial temporal storage
    """
    print('Simulation...', end=''); tim = time.time()

    ### THIS IS THE PART WHERE YOU CAN PUT SOME COMMUNICATION
    # :::C::: LAUNCH THE OTHER PROGRAM

    for i in range(param['Nt']-1):
        y += M.rk4(y, op, param)                       # The vector y is dynamically updated 
        t += param['dt']
        tprevious += param['dt']          # time is incremented
        if tprevious >= param['Tstore']:              # if it is time to record the state
            idx+=1                                # we give it a new index
            Y_s[:, idx, :] = np.copy(y)             # we write it in the "book" Y_s
            t_s[idx]     = t*1                    # we write the time 
            tprevious=0.                          # we reinitialise time before next record

        ### THIS IS THE PART WHERE YOU CAN PUT SOME COMMUNICATION
        """
        :::C::: The elements for communications between programs.

        The typical method is to put a value (or a list) in a .txt file, make it read by the other file

        The important thing: You have to wait one program to finish and write the file for the other one to start
        Two rough strategies:
            * On each side check the date of last update for each communication 
            * Have a common writing flag with a flag (if value is "A" then it's python, "B" is vensim)

        The elegant strategies:
            * Using a python-Vensim binding with step-by-step Dymends in a big function 
            * Using python to launch some Shell-type code

        f = open("demofile2.txt", "w")
        f.write(Content)
        f.close()

        f = open("demofile2.txt", "r")
        Content = f.read()
        f.close()

        import os.path, time
        print("Last modified: %s" % time.ctime(os.path.getmtime("test.txt")))
        print("Created: %s" % time.ctime(os.path.getctime("test.txt")))

        For a better version, this has to be moved into the dynamic function f_
        """

    print('done ! elapsed time :', time.time()-tim,'s')

    # Save the data as a pickle file in a new folder
    if param['Save']:
        FG.savedata(rootfold, t, Y_s, param, op)

    # ----------------------
    # Results interpretation
    """
    Now that the simulation is done, we can translate its results in a more readable fashion.
    r is the expansion of Y_s into all the relevant variables we are looking for, stored as a dictionnary.
    Then, the other parts are simply plots of the result
    """

    r = FG.expandY(Y_s, t_s, op, param)  # Result dictionnary 
    r = M.getperiods(r, param, op)      # Period measurements 

    # GRAPHS, 
    if plot is True:
        plts.PhilAndInvest(param,)             # Show behavior functions 
        plts.GoodwinKeenTypical(r, param,)      # Typical 3-Dimension phase-plot
        plts.omegalambdacycles(r, param,)       # 2-D omega-lambda phase portrait
        plts.GraphesIntensive(r, param)
        plts.GraphesExtensive(r, param)
        #plts.PhasewithG(r, param, op,)           # Phase diagram with growth
        plts.PeriodPlots(r, param, op,)          # Study stability-period
        #plts.map2DLambdaT(r, op,)
        #plts.MesoMeanSTD(r, param, op,)

    return



# #############################################################################
# #############################################################################
#                   Handle bash
# #############################################################################


if __name__ ==  '__main__':

    # executable help
    msg = """Thus function does this and that"""

    # Instanciate parser                                                        
    parser = argparse.ArgumentParser(description=msg)

    # Define input arguments                                                                                 
    parser.add_argument(
        '-p',
        '--plot',
        type=_str2bool,
        help='flag indicating whether to plot figures',
        default=_PLOT,
        required=False,
    )
    parser.add_argument(
        '-root',
        '--rootfold',
        type=str,
        help='path where to save data',
        default=_ROOTFOLD,
        required=False,
    )
    kwdargs = dict(parser.parse_args()._get_kwargs())

    # Call function                                                             
    run(**kwdargs)
