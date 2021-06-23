# -*- coding: utf-8 -*-
"""
PARAMETERS FOR GOODWIN-KEEN LIKE MODELS
"""
 
import numpy as np
import FunctionsGoodwin as FG

def initparams():
    """
    Create a dictionnary containing all the parameters necessary for the simulation. 
    Their description is in comment here. 
    """
    
    
    p = {} # Dictionnary of parameters 
    # SET OF EQUATIONS ########################################################
    p['EquationSet'] =FG.f_reduced           # The system of equation the code will solve
    
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
    
    ### Numerical PARAMETERS ##################################################
    # Time Discretisation ###
    p['Tmax']       = 100                   # Time simulated (in years)
    p['time_safety']= 30                    # Factor of time increment reduction. The system calculate an increment dtmax, then divide it by this value
    
    # Space discretisation/Number of systems 
    p['Nx']  = 2                           # Number of parrallel economy (that can be coupled)
    
    # DATA STORAGE #####
    p['StorageMode']='Full'                 # Keep all iterations of only few
    p['Tstore']     = 0.1                   # Time between two storage (if StorageMode=full, it goes to dt)
    p['Save']       = False 
    p['printblabla']= True
    
    
    
    ### EACH MODULES AND VALUES ###############################################
    # Population evolution 
    p['beta']   = 0.025                      # Rate of population growth     (y^-1)
    p['PopSat'] = 10000                      # Ratio between initial and maximal population 
    p['alpha' ] = 0.02                       # Rate of productivity increase (y^-1)
    
    # Capital properties
    p['delta2'] = 0.005                      # Rate of ALL capital depletion     (y^-1)
    p['delta1'] = 0.005                      # Rate of WORKING capital depletion (y^-1)
    p['pK']     = 0                          # Price consumption of capital
    
    p['nu'   ] = 3                           # Kapital to output ratio as in Leontiev. CAREFUL IN CES this is 1/A
    p['eta']   = 10000                       # CES Only 1/(1+substituability)
    p['b']     = .5                          # CES parameter : capital part of the production 
    p['z']     = 1                           # Markup on salary estimation
    
    # INTEREST / Price
    p['r']        = .03  # Interest at the bank
    p['eta']     = .192 # Typical rate for inflation
    p['mu']       = 1.3  # Mark-up of price
    p['gamma']   = 1 # Money-illusion 
    
    
    
    ### FUNCTIONS AND THEIR PARAMETERS ########################################
    # PHILIPS CURVE (employement-salary increase)
    p['phiType']  = 'divergent'              # 'divergent' , 'linear' modes
    p['phinul']   = 0.04                     # Unemployement rate at which there is no salary increase with no inflation
    p['phislope'] = 2                       # slope in a linear model (power in a diverging model)         
    
    # KEEN INVESTMENT FUNCTION (profit-investment function)
    p['k0'] = -0.0065
    p['k1'] = np.exp(-5)
    p['k2'] = 20
    
    # LINEAR DIVIDENT PROFITS 
    p['div0'] = 0.138 # Part of GDP as dividends when pi=0
    p['div1'] = 0.473 # Slope 
    
    return p


    """
    ################################################################################################
    # Coupling Effets ##############################################################################
    p['g1']  = .0                          # GLOBAL         EFFECTS OF LAMBDA
    p['g2']  = .00                          # WITH NEIGHBORS EFFECTS OF LAMBDA
    p['muI'] = 0                            # 
    p['muN'] = 0                            #
     
    
    # Recrutment-Firing dynamics
    #tau = np.linspace(0.8,1.1,p['Nx'])
    p['tauR'] = 0.0#85                      # Typical time for recruitement (y)
    p['tauF'] = 0.0                         # Typical time for firing       (y)
    
    p['tauL']   = 0.0                       # Typical relaxation time on employement perception
    p['tauLam'] = 2.
    p['tauK']   = 2                   # Typical relaxation time on  
        
    
        
    
    ################################################################################################
    # CLIMATE-INTEGRATION ##########################################################################
    '''
    Parameters from Debt&Damages
    '''
        
    p['theta'] = 2.6                     # Convexity on abattement cost function
    p['dsigma'] = - 0.001               # Variation rate of the growth of emission intensity 
    p['dPBS']   = - 0.005               # Growth rate of back-stop technology price   
    p['dPc']   = 0 #[CORRECT]               # Growth rate of back-stop technology price 
    p['dEland'] = - 0.022               # Growth rate of land use change in CO2 emission
    
    # Damage function (on GDP)
    '''
    D = 1 - (1 + p['pi1']*T + p['pi2']*T**2 + p['pi3']*T**p['zeta'] )**(-1)
    '''
    p['pi1'] = 0                         # Linear temperature impact
    p['pi2'] = .00236                    # Quadratic temperature impact
    p['pi3'] = .00000507                 # Weitzmann Damage temperature impact 
    p['zeta']= 6.754                     # Weitzmann impact 
    p['fk']  = 1/3                       # Fraction of environmental damage allocated to the stock of capital
    
    # Climate model
    p['Phi12']   = .024 #Transfer of carbon from atmosphere to biosphere
    p['Phi23']   = .001 #Transfer from biosphere to stock
    
    p['C' ]    = 1/.098 # Heat capacity of fast-paced climate
    p['C0']    = 3.52   # Heat capacity of inertial component of climate
    p['gamma'] = 0.0176 # Heat exchange coefficient between layer
    p['Tsens'] = 3.1    # Climate sensibility
    
    p['dFexo'] = 0
    p['FexoMax'] = 0 
    p['F2CO2']   = 3.681 # W/M2, doubling CO2 impact on forced radiations
       
    #####################################################################
    ### INITIAL CONIDITONS 
    p['CO2at_ini'] = 851
    p['CO2up_ini'] = 460
    p['CO2lo_ini'] = 1740
    p['T_ini'] = 0
    """
    


