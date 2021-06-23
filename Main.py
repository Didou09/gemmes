#!/usr/bin/env python
# -*- coding: utf-8 -*-


#rootfold='/home/ejp/Desktop/Goodwin-Keen/simulations'
#rootfold='D:\\Georgetown\\Simulation-Goodwin\\'
rootfold = os.path.dirname(__file__)

### WELCOME MESSAGE ################################################################################
"""
GOODWIN-TYPE RESOLUTION ALGORITHM 

This python code has been created to simulate Goodwin-Keen simulations, with all the variants to it. 
Do not worry about the number of section : the architecture is thought to be user-friendly and coder-friendly

WHAT ARE THE SPECIFICITIES OF THIS CODE ? 
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

HOW TO IMPLEMENT YOUR TOY-MODEL ? 
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

### INITIALISATION LIBRAIRIES ######################################################################
####################################################################################################
for _ in range(1):
    """
    The code is importing some other libraries so that we can use their functions. 
    """
    import numpy as np              # Array/Numerical operations
    import time                     # Time (run speed) printing
    import matplotlib.pyplot as plt # Plot elements 
    plt.close('all')                # Remove to 
    
    """
    These libraries are local code fragments we use later
    """
    import Parameters as Par      # All the parameters of the system
    import FunctionsGoodwin as FG # Core of models
    import plots as plts          # Already written plot functions
    import Miscfunc as M          # All miscellaneous functions

    ### MATPLOTLIB INTERACTIVITY ####
    '''
    Matplotlib allow in-line (in terminal) output or interactive output. Select one here.
    You can change it in the jupyter (interactive terminal) of your IDE also 
    %matplotlib # Allow interactive output
    %matplotlib inline # plot in terminal
    '''
    #%matplotlib inline 

print('List of existing set of equations  :')
for f in [ f for f in dir(FG) if 'f_' in f] : print(f)
print(3*'\n')

### PARAMETERS INITIALISATION ######################################################################
####################################################################################################    
for _ in range(1):
    """
    This part create 'p' (parameters dic) and 'op' (operators dic) depending of the system size and properties
    Typically, p contains the "default parameters" values and should be open in another windows, 
    I edit them after p initialisation if I want to do some original runs 
    """
    
    p = Par.initparams()              # PARAMETERS IN THIS FILE ( PARAMETERS )
    
    p['lambdamax']  = 1-10**(-2)          # Variation of lambda on the simulation 
    
    # The good place to add a loop on certain parameters if you need 
    p = M.preparePhilips(p)    # Determining Philips curve and Solow point
    p = M.preparedT(p)         # Determine the best value for dt and associates
    op= []#M.prepareOperators(p)  # Spatial operators initialisation

### initial conditions #############################################################################
####################################################################################################
"""
The initial state vector at t=0 of the system in a "readable" language.
Depending of the equation set you use, all of them will not necessary be used. 
(example, a reduced set will read ic['d'], a complete set will read ic['D']) 
"""

ic={} # Dictionnary of initial conditions
v1 = np.ones(p['Nx']) ### Pratical notation for more readable code : an array with 1 everywhere

ic['a']       = 1.*v1            # Initial productivity per worker
ic['N']       = 1.*v1            # Initial quantity of people

ic['d']       = 0*1.53*v1 
ic['omega']   = 1.*v1*p['omega0']#*np.arange(0,1,p['Nx'])
ic['lambda']  = np.linspace(p['lambdamin'],p['lambdamax'],p['Nx'])#np.linspace(p['lambdamin'],p['lambdamax'],p['Nx'])    # Distribution of initial lambda
ic['L']       = ic['omega']*ic['N']
ic['K']       = ic['a']*ic['L']*p['nu']                                # Beginning at equilibrium
ic['W']       = p['omega0'] *ic['K']/p['nu']                           # Initial distribution of salary


'''
ic['D']     =1.*v1            # Debt   
ic['sigma'] =1.*v1            # Intensity of emission in production
ic['P']     =1.*v1            # Global price
ic['Pbs']   =547.22*v1            # Price Backstop
ic['Pc']    =1.*v1            # Carbon price


ic['Fexo']  =1.*v1            # Exogenous forcing intensity
ic['Eland'] =2.6*v1           # Emission from land    
ic['CO2at'] =1.*v1*p['CO2at_ini']           # CO2 in atmosphere
ic['CO2up'] =1.*v1*p['CO2up_ini']           # CO2 in upper ocean and biosphere
ic['CO2lo'] =1.*v1*p['CO2lo_ini']          # CO2 in deep ocean
ic['T']     =1.*v1           # Atmospheric temperature 
ic['T0']    =1.*v1           # Deeper ocean temperature 
'''
ic['t']     =1.*v1           # Year in the simulation 


y,p = FG.prepareY(ic,p) ### The vector y containing all the state of a time t.
### Initialisation of data storage #################################################################
####################################################################################################
"""
This part is very close to the numerical resolution core. No need to delve in it
It looks "complex" because it allows the record of few timestep only
"""
Y_s = np.zeros((p['Nvar'],p['Ns'],p['Nx']))     # stock dynamics
Y_s[:,0,:]= 1*y                                 # first stock
t         = 0
t_s       = np.zeros(p['Ns'])                   # stock time
tprevious = 0                                   # deltatime between two records
idx       = 0                                   # index of the current record

### Simulation iteration ############################################################################
#####################################################################################################
"""
The RK4 core, with partial temporal storage 
"""
print('Simulation...',end='');tim=time.time()

### THIS IS THE PART WHERE YOU CAN PUT SOME COMMUNICATION
# :::C::: LAUNCH THE OTHER PROGRAM

for i in range(p['Nt']-1):
    y += M.rk4(y,op,p)                       # The vector y is dynamically updated 
    t += p['dt'];tprevious+= p['dt']          # time is incremented
    if tprevious >= p['Tstore']:              # if it is time to record the state
        idx+=1                                # we give it a new index
        Y_s[:,idx,:] = np.copy(y)             # we write it in the "book" Y_s
        t_s[idx]     = t*1                    # we write the time 
        tprevious=0.                          # we reinitialise time before next record
        
    ### THIS IS THE PART WHERE YOU CAN PUT SOME COMMUNICATION
    """
    :::C::: The elements for communications between programs.
    
    The typical method is to put a value (or a list) in a .txt file, make it read by the other file
    
    The important thing : You have to wait one program to finish and write the file for the other one to start
    Two rough strategies : 
        * On each side check the date of last update for each communication 
        * Have a common writing flag with a flag (if value is "A" then it's python, "B" is vensim)
    
    The elegant strategies : 
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

if p['Save'] : FG.savedata(rootfold,t,Y_s,p,op) # Save the data as a pickle file in a new folder
### Results interpretation #########################################################################
####################################################################################################
"""
Now that the simulation is done, we can translate its results in a more readable fashion. 
r is the expansion of Y_s into all the relevant variables we are looking for, stored as a dictionnary.
Then, the other parts are simply plots of the result
"""

r = FG.expandY(Y_s,t_s,op,p)  # Result dictionnary 
r = M.getperiods(r,p,op)      # Period measurements 

# GRAPHS, 
plts.PhilAndInvest(p,)             # Show behavior functions 
plts.GoodwinKeenTypical(r,p,)      # Typical 3-Dimension phase-plot
plts.omegalambdacycles(r,p,)       # 2-D omega-lambda phase portrait
plts.GraphesIntensive(r,p)    
plts.GraphesExtensive(r,p) 
#plts.PhasewithG(r,p,op,)           # Phase diagram with growth
plts.PeriodPlots(r,p,op,)          # Study stability-period
#plts.map2DLambdaT(r,op,)
#plts.MesoMeanSTD(r,p,op,)
# :-)
# :-)
# :-)
# :-)
# :-)
# :-)
