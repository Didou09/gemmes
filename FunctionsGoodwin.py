import numpy as np
import Miscfunc as M

###############################################################################
### CORE RESOLUTION ###########################################################
###############################################################################
def f_reduced(y,op,p):
    """
    Keen95 system
    + money illusion with inflation
    """
    omeg = y[0]
    lamb = y[1]
    d    = y[2]
    
    pi = 1 - omeg - p['r']*d # pi*Y is the net profit of the enterprise
    i = p['eta']*(omeg*p['mu']-1)
    
    
    op = omeg * (M.philips(lamb,p)-p['alpha']-p['gamma']*i)
    lp = lamb * (M.fk(pi,p)/p['nu'] - p['alpha'] - p['alpha'] - p['beta'] )
    dp = M.fk(pi,p) - pi - d*  (M.fk(pi,p)/p['nu'] - p['delta1'] +i)
    return np.array([op*p['dt'],lp*p['dt'],dp*p['dt']])

def f_ces_reduced(y,op,p):
    """
    Set of equations for Goodwin-Keen,a CES, and salary defined as marginal return on labour . 
    From : https://doi.org/10.1007/s11579-018-0231-6
    
    A from the article is here substitued by 1/nu from Leontiev"""
    omeg = y[0]
    lamb = y[1]
    d    = y[2]

    pi = 1 - omeg - p['r']*d
    op = omeg * (p['eta']/(1+p['eta']))*(M.philips(lamb,p)-p['alpha'])   
    lp = lamb * ( M.fk(pi,p)/p['nu']* ((1-omeg )/p['b'])**(1/p['eta'])- op / (p['eta']*(1-omeg)*omeg) -p['delta1']-p['alpha']-p['beta'] )
    dp = d    * ( p['r'] - M.fk(pi,p)/p['nu']* ((1-omeg )/p['b'])**(1/p['eta']) + p['delta1'] +  op / (p['eta']*(1-omeg))) + M.fk(pi,p)-(1-omeg)
    return np.array([op*p['dt'],lp*p['dt'],dp*p['dt']])

def f_GEMMES_CLIM(y,op,p):
    """
    This system only use the climate module of GEMMES, this is a test to check everything is working
    """
    Fexo  = y[0]           # Exogenous forcing intensity
    Eland = y[1]           # Emission from land    
    CO2at = y[2]           # CO2 in atmosphere
    CO2up = y[3]           # CO2 in upper ocean and biosphere
    CO2lo = y[4]           # CO2 in deep ocean
    T     = y[5]           # Atmospheric temperature 
    T0    = y[6]           # Deeper ocean temperature 
 
    F = Fexo 
    E = Eland
       
    CO2atup = CO2at / CO2up  # Ratio of CO2 concentration (chemical potential effect ?)
    CO2uplo = CO2up / CO2lo  # Ratio of CO2 concentration (chemical potential effect ?)
    
 
    CO2atp = E + CO2at*(-p['Phi12']) + CO2up*( p['Phi12']*CO2atup)            + CO2lo * 0
    CO2upp = 0 + CO2at*( p['Phi12']) + CO2up*(-p['Phi12']*CO2atup-p['Phi23']) + CO2lo * (+p['Phi23']*CO2uplo) 
    CO2lop = 0 + CO2at*( 0         ) + CO2up*( 0                 +p['Phi23']) + CO2lo * (-p['Phi23']*CO2uplo)  
     
    Elandp = Eland * p['dEland']    
    Fexop  = Fexo * p['dFexo']*(1-Fexo/p['FexoMax'])    

    Tp  = (F-p['rho']*T - p['gamma']*(T-T0) )/p['C'] # Atmosphere temperature evolution : forcing, dissipation in space, transfer to lower strates
    T0p = p['gamma']*(T-T0)/p['C0']                  # Lower strates temperature evolution 

      
    return np.array([
                    Fexop  *p['dt'],           # Exogenous forcing intensity
                    Elandp *p['dt'],           # Emission from land    
                    CO2atp *p['dt'],           # CO2 in atmosphere
                    CO2upp *p['dt'],           # CO2 in upper ocean and biosphere
                    CO2lop *p['dt'],           # CO2 in deep ocean
                    Tp     *p['dt'],           # Atmospheric temperature 
                    T0p    *p['dt'],           # Deeper ocean temperature 
                            p['dt']])          # Time evolution

def f_GEMMES_CORE(y,op,p):
    """
    Just the economic module with no climate retroaction 
    """
    K     = y[0]           # Capital 
    a     = y[1]           # Productivity (written A elsewhere)
    N     = y[2]           # World population
    W     = y[3]           # Salary 
    P     = y[4]           # Global price
    D     = y[5]           # Debt

    T     = 0
    sigma = y[6]           # Intensity of emission in production
    Pbs   = y[7]           # Price Backstop
    Pc    = y[8]           # Carbon price    
    t     = y[-1]
   ### 2) We determine the accounting quantities ################################################## 
    ## 2.1) Economy 
    # Damage     
    D  =  1 - (1 + p['pi1']*T + p['pi2']*T**2 + p['pi3']*T**p['zeta'])**(-1) # Damage function from Dietz and Stern
    Dk = p['fk']*D
    Dy = 1 - (1-D)/(1-Dk)
    deltak = p['delta1']+Dk
    
    # Production
    Y0 = K / p['nu']                                           # GDP without removal (Leontiev optimized)
    L = K /(a*p['nu'])                                         # Equivalent labor (Leontiev optimized)

    #n = np.minimum( (Pc/Pbs  )**(1/(p['theta']-1)), 1) # Edogeneous emission reduction rate 
    n=1
    A  = n**(p['theta'])/p['theta'] * sigma*Pbs                # Abbatement 
    Eind = Y0*sigma*(1-n)                                      # Emission induced by the economy
    Y = Y0 * ( 1 - Dy)*(1-A)                                   # REAL GDP

    # Paying carbone
    Eind = Y0*sigma*(1-n)                                      # Emission induced by the economy
    Ct = Pc*Eind                                               # Total in Carbon tax unit 

    # Production repartition 
    PI =(P*Y          + # Product sold
        (-W*L)        + # Wage paid 
        (-p['r']*D)   + # Interest paid
        (-P*K*deltak) + # THERE IS A DEPRECIATION OF CAPITAL HERE ! 
        (-P*Ct))        # Carbon Tax payment
    
    pi  = PI/(P*Y)           # Profit ratio of productin
    PId = M.delta(pi,p)*P*Y  # Shareholder Profit 
    PIr = PI - PId           # Real profit

    if p['CostUnit']=='DebtDamage': c =  W*L/(P*Y)                # Cost per unit of production. VERSION DEBT&DAMAGE
    elif p['CostUnit']=='Coping' :  c = (((1-A)(1-Dy) )**(-1) *   # VERSION COPING
                                        (W/(P*a) + p['nu'] * (p['r']*D+PId)/(p*K) + Pc*sigma*(1-n)))

    i = p['etaP']*(P*p['mu']*c-1) # Price inflation
    I = Y*M.fk(pi,p)                # Investment
     
    ### CALCULATION OF SYSTEM EVOLUTION IN TIME ####################################################
    ### 3.1.a) We calculate variations (economy)    
    ap = a * p['alpha']                                 # Evolution of labour efficiency (exogenous technical progress)
    Np = N * p['beta']  * (1- N/p['PopSat'])            # Evolution of population (Logistic function)   
    Pp = P * i                                          # Inflation (Unstable Relaxation)
    Wp = W * M.philips(L/N,p)   -p['gammaP']*i          # Salary evolution (Instant)
    Kp =-K * deltak + I                                 # Evolution of capital    
    Dp = (P*I +                                         #
          - PIr)
   
    # 3.1.b) Economy-climate price and intensity 
    sigmap = p['dsigma'] * sigma    
    Pbsp   = p['dPBS']*Pbs                               # Diminution of the price for backstop technology
    Pcp    = p['dPc'] *Pc                               # Purely exponential dynamics
    #Pcp    = Pc * (p['aPc']+p['bPc']/t)                  # With a short-term boost         

    return np.array([
                    Kp     *p['dt'],           # Capital 
                    ap     *p['dt'],           # Productivity
                    Np     *p['dt'],           # World population
                    Wp     *p['dt'],           # Number of workers
                    Pp     *p['dt'],           # Global price
                    Dp     *p['dt'],           # Debt 
                    sigmap *p['dt'],           # Intensity of emission in production
                    Pbsp   *p['dt'],           # Price Backstop
                    Pcp    *p['dt'],           # Carbon price
                            p['dt']])          # Time evolution
          
def f_GEMMES(y,op,p):
    """
    The equations extracted from "Debt & Damages"
    https://doi.org/10.1016/j.inteco.2018.02.002
    """
    
    ### 1) We transform y (the vector state) into variable names ###################################
    
    # 1.a) Pure economy
    K     = y[0]           # Capital 
    a     = y[1]           # Productivity (written A elsewhere)
    N     = y[2]           # World population
    W     = y[3]           # Salary 
    P     = y[4]           # Global price
    D     = y[5]           # Debt
    
    # 1.b) Economy related to emission
    sigma = y[6]           # Intensity of emission in production
    Pbs   = y[7]           # Price Backstop
    Pc    = y[8]           # Carbon price

    # 1.c) CO2 and temperature module
    Fexo  = y[-8]           # Exogenous forcing intensity
    Eland = y[-7]           # Emission from land    
    CO2at = y[-6]           # CO2 in atmosphere
    CO2up = y[-5]           # CO2 in upper ocean and biosphere
    CO2lo = y[-4]           # CO2 in deep ocean
    T     = y[-3]           # Atmospheric temperature 
    T0    = y[-2]           # Deeper ocean temperature 
    
    t     = y[-1]           # Time 
    
    ### 2) We determine the accounting quantities ################################################## 
    ## 2.1) Economy 
    # Damage 
    D  =  1 - (1 + p['pi1']*T + p['pi2']*T**2 + p['pi3']*T**p['zeta'])**(-1) # Damage function from Dietz and Stern
    Dk = p['fk']*D
    Dy = 1 - (1-D)/(1-Dk)
    deltak = p['delta1']+Dk
    
    # Production
    Y0 = K / p['nu']                                           # GDP without removal (Leontiev optimized)
    L = K /(a*p['nu'])                                         # Equivalent labor (Leontiev optimized)

    n = np.minimum( (Pc/Pbs  )**(1/(p['theta']-1)), 1 ,axis=0) # Edogeneous emission reduction rate 
    A  = n**(p['theta'])/p['theta'] * sigma*Pbs                # Abbatement 
    Eind = Y0*sigma*(1-n)                                      # Emission induced by the economy
    Y = Y0 * ( 1 - Dy)*(1-A)                                   # REAL GDP

    # Paying carbone
    Eind = Y0*sigma*(1-n)                                      # Emission induced by the economy
    Ct = Pc*Eind                                               # Total in Carbon tax unit 

    # Production repartition 
    PI =(p*Y          + # Product sold
        (-W*L)        + # Wage paid 
        (-p['r']*D)   + # Interest paid
        (-p*K*deltak) + # THERE IS A DEPRECIATION OF CAPITAL HERE ! 
        (-p*Ct))        # Carbon Tax payment
    
    pi  = PI/(P*Y)           # Profit ratio of productin
    PId = M.delta(pi,p)*p*Y  # Shareholder Profit 
    PIr = PI - PId           # Real profit

    if p['CostUnit']=='DebtDamage': c =  W*L/(P*Y)                # Cost per unit of production. VERSION DEBT&DAMAGE
    elif p['CostUnit']=='Coping' :  c = (((1-A)(1-Dy) )**(-1) *   # VERSION COPING
                                        (W/(P*a) + p['nu'] * (p['r']*D+PId)/(p*K) + Pc*sigma*(1-n)))

    i = p['etaP']*(P*p['mu']*c-1) # Price inflation
    I = Y*M.fk(pi,p)                # Investment
    
    ### 2.2 Climate Module
    Eland = Eland            # Emission induced by nature. WRONG FORMULA 
    E = Eind + Eland         # Total emissions

    CO2atup = CO2at / CO2up  # Ratio of CO2 concentration (chemical potential effect ?)
    CO2uplo = CO2up / CO2lo  # Ratio of CO2 concentration (chemical potential effect ?)
    
    Find = (p['F2CO2'] /np.log(2)) * np.log( CO2at / p['CO2at_ini'] ) # Forcing created by industry
    F = Find + Fexo
    
    ### CALCULATION OF SYSTEM EVOLUTION IN TIME ####################################################
    ### 3.1.a) We calculate variations (economy)    
    ap = a * p['alpha']                                 # Evolution of labour efficiency (exogenous technical progress)
    Np = N * p['beta']  * (1- N/p['PopSat'])            # Evolution of population (Logistic function)   
    Pp = P * i                                          # Inflation (Unstable Relaxation)
    Wp = W * M.philips(L/N,p)   -p['gammaP']*i          # Salary evolution (Instant)
    Kp =-K * deltak + I                                 # Evolution of capital    
    Dp = (p*I +                                         #
          - PIr)
   
    # 3.1.b) Economy-climate price and intensity 
    sigmap = p['dsigma'] * sigma    
    Pbsp   = p['dPBS']*Pbs                               # Diminution of the price for backstop technology
    #Pcp    = p['dPc'] *Pc                               # Purely exponential dynamics
    Pcp    = Pc * (p['aPc']+p['bPc']/t)                  # With a short-term boost


    # 3.2  CO2 3-Layers Dynamics 
    CO2atp = E + CO2at*(-p['Phi12']) + CO2up*( p['Phi12']*CO2atup)            + CO2lo * 0
    CO2upp = 0 + CO2at*( p['Phi12']) + CO2up*(-p['Phi12']*CO2atup-p['Phi23']) + CO2lo * (+p['Phi23']*CO2uplo) 
    CO2lop = 0 + CO2at*( 0         ) + CO2up*( 0                 +p['Phi23']) + CO2lo * (-p['Phi23']*CO2uplo)  
     
    Elandp = Eland * p['dEland']    
    Fexop  = Fexo * p['dFexo']*(1-Fexo/p['FexoMax'])    
    
    Tp  = (F-p['F2CO2']*T/p['Tsens'] - p['gamma']*(T-T0) )/p['C'] # Atmosphere temperature evolution : forcing, dissipation in space, transfer to lower strates
    T0p = p['gamma']*(T-T0)/p['C0']                  # Lower strates temperature evolution 


    ### 4) The variation goes back into the code ! #################################################   
    return np.array([
                    Kp     *p['dt'],           # Capital 
                    ap     *p['dt'],           # Productivity
                    Np     *p['dt'],           # World population
                    Wp     *p['dt'],           # Number of workers
                    Pp     *p['dt'],           # Global price
                    Dp     *p['dt'],           # Debt 
                    sigmap *p['dt'],           # Intensity of emission in production
                    Pbsp   *p['dt'],           # Price Backstop
                    Pcp    *p['dt'],           # Carbon price
                    Fexop  *p['dt'],           # Exogenous forcing intensity
                    Elandp *p['dt'],           # Emission from land    
                    CO2atp *p['dt'],           # CO2 in atmosphere
                    CO2upp *p['dt'],           # CO2 in upper ocean and biosphere
                    CO2lop *p['dt'],           # CO2 in deep ocean
                    Tp     *p['dt'],           # Atmospheric temperature 
                    T0p    *p['dt'],           # Deeper ocean temperature 
                            p['dt']])          # Time evolution
                   
def f_reduced_relax(y,op,p):
    """
    Goodwin-Keen system with relaxation on capital, on employement. Two depreciation terms depending on K0 or K
    Equation : Paul report march 10
    """    
    g           = y[0]
    lambdapoint = y[1]
    lambd       = y[2]
    omega       = y[3]
    d           = y[4]
    
    pi = 1 - omega - p['r']*d
    gp  = (1/p['tauK'])*(M.fk(pi,p)/p['nu'] - p['delta1']*(1-p['tauK']*p['delta2']) - g*(1+p['tauK']*(p['delta1']+p['delta2'])))
    lpp = (lambd/p['tauLam'])*(g-p['beta']-p['alpha']-lambdapoint)
    lp  = lambdapoint
    op  = omega*(M.philips(lambd,p)-p['alpha'])
    dp  = M.fk(pi,p)-pi- g
    
    return np.array([gp *p['dt'] , 
                     lpp*p['dt'] , 
                     lp *p['dt'] , 
                     op *p['dt'] , 
                     dp *p['dt'] ])    
    
def f_full_leontiev(y,op,p):
    '''
    ALL THE HYPOTHESIS, IT'S A BIG SYSTEM
    We keep all equations explicit so that the system is easier to tweak
    1) We transform the vector state into the name of variables (easier to read)
    2) We calculate intermediary quantity
    3) We deduce the time variation for differential equations
    4) We send the variation through an instant dt
    '''
    
    ### 1) We transform y (the vector state) into variable names
    a = y[0]
    K = y[1]
    L = y[2]
    W = y[3]
    N = y[4]
    D = y[5]
    
    # LEONTIEV
    Y = K / p['nu']                 # Production function
    L0 = K / (a*p['nu'])             # Labour instant adaptation

    ### 2.1 we calculate the type of relaxation on L 
    L=L0
    Lp = (L0-L)/p['dt']
    '''signL = np.sign(L-L0)
    instant=( (signL== 1)*(p['tauF']==0) + 
              (signL==-1)*(p['tauR']==0) + 
              (signL== 0)).astype(np.bool) # True if L is at equilibrium    
    Tautemp = (1-instant)* ( (signL+1)/2 * p['tauF'] + 
                             (1-signL)/2 * p['tauR'] )+ \
              (  instant)* (      1    ) * p['dt']
    L[instant]=L0[instant]
    Tautemp[instant]=p['dt']
    Lp = (L0-L)/Tautemp    
    '''
    ### 2.2 We apply it to the system
    P = Y - L*W - p['r']*D#- p['pK']*K  # Profit
    I = M.fk(P/Y,p)*Y                 # Investment
    lambd_local  = L/N
    #lambd_meso   = (L + np.roll(L,1) + np.roll(L,-1))/(N + np.roll(N,1) + np.roll(N,-1))
    #lambd_global = np.dot(op['Mean'],L)/np.dot(op['Mean'],N)

    ### 3) We calculate variations    
    Kp =-K * p['delta1'] + I         # Evolution of capital
    Ap = a * p['beta']              # Evolution of labour efficiency
    Np = N * p['alpha']             # Evolution of population
    Dp = I - P                      # Debt evolution
    Wp = W * M.philips(L/N,p)         # Salary evolution
    

    #Lp = (L0-L)/Tautemp          
    Wp =W*((1-p['g1']-p['g2'])*M.philips( lambd_local ,p))
          #+(  p['g1']        )*philips( lambd_global,p)
          #+(          p['g2'])*philips( lambd_meso  ,p))

    ### 4) The variation goes back into the code !    
    return np.array([Ap*p['dt'] , 
                     Kp*p['dt'] , 
                     Lp*p['dt'] , 
                     Wp*p['dt'] , 
                     Np*p['dt'] ,
                     Dp*p['dt'] ])    
    
def f_full_CES(y,op,p):
    '''
    ALL THE HYPOTHESIS, IT'S A BIG SYSTEM
    We keep all equations explicit so that the system is easier to tweak
    1) We transform the vector state into the name of variables (easier to read)
    2) We calculate intermediary quantity
    3) We deduce the time variation for differential equations
    4) We send the variation through an instant dt
    '''
    
    ### 1) We transform y (the vector state) into variable names
    a = y[0]
    K = y[1]
    L = y[2]
    W = y[3]
    N = y[4]
    D = y[5]
    
    ### 2) we calculate intermediary quantities        
    # CES 
    Yc = (1/p['nu'])*K*p['b']**(-1/p['eta'])           # Characteristic value of capital
    Lc = K/a * ((1-p['b'])/p['b'])**(1/p['eta'])       # Characteristic quantity of labor
    Y = Yc * (1 + (L/Lc)**(-p['eta']))**(-1/p['eta'])  # GDP
    
   

    ### 2.1 we calculate the type of relaxation on L 
    L0 = Lc * (p['z']*Lc*W/Yc)**(-1/(p['eta']+1))      # BE CAREFUL, HERE WE ARE NOT TAKING THE EXACT VALUE.
 
    signL = np.sign(L-L0)
    instant=( (signL== 1)*(p['tauF']==0) + 
              (signL==-1)*(p['tauR']==0) + 
              (signL== 0)).astype(np.bool) # True if L is at equilibrium    
    Tautemp = (1-instant)* ( (signL+1)/2 * p['tauF'] + 
                             (1-signL)/2 * p['tauR'] )+ \
              (  instant)* (      1    ) * p['dt']
    L[instant]=L0[instant]
    Tautemp[instant]=p['dt']
    Lp = (L0-L)/Tautemp    
    
    ### 2.2 We apply it to the system
    P = Y - L*W - p['r']*D#- p['pK']*K  # Profit
    I = M.fk(P/Y,p)*Y                 # Investment
    lambd_local  = L/N
    #lambd_meso   = (L + np.roll(L,1) + np.roll(L,-1))/(N + np.roll(N,1) + np.roll(N,-1))
    #lambd_global = np.dot(op['Mean'],L)/np.dot(op['Mean'],N)

    ### 3) We calculate variations    
    Kp =-K * p['delta1'] + I         # Evolution of capital
    Ap = a * p['beta']              # Evolution of labour efficiency
    Np = N * p['alpha']             # Evolution of population
    Dp = I - P                      # Debt evolution
    Wp = W * M.philips(L/N,p)         # Salary evolution
    

    #Lp = (L0-L)/Tautemp          
    Wp =W*((1-p['g1']-p['g2'])*M.philips( lambd_local ,p))
          #+(  p['g1']        )*philips( lambd_global,p)
          #+(          p['g2'])*philips( lambd_meso  ,p))

    ### 4) The variation goes back into the code !    
    return np.array([Ap*p['dt'] , 
                     Kp*p['dt'] , 
                     Lp*p['dt'] , 
                     Wp*p['dt'] , 
                     Np*p['dt'] ,
                     Dp*p['dt'] ])    
   
##############################################################################
### SYSTEM INITIALISATION ###############################################
###############################################################################
def prepareY(ic,p):
    """
    Transform the dictionnary ic into an array y
    """
    if (p['EquationSet']==f_reduced) or (p['EquationSet']==f_ces_reduced):
        p['Nvar'] = 3
        y = np.zeros((p['Nvar'],p['Nx']))
        
        y[0,:] = ic['omega']  
        y[1,:] = ic['lambda'] 
        y[2,:] = ic['d']      
                
    elif (p['EquationSet']==f_full_leontiev) or (p['EquationSet']==f_full_CES):
        p['Nvar'] = 6
        y = np.zeros((p['Nvar'],p['Nx']))
        
        y[0,:] = ic['a']  
        y[1,:] = ic['K'] 
        y[2,:] = ic['L'] 
        y[3,:] = ic['W'] 
        y[4,:] = ic['N'] 
        y[5,:] = ic['D']          
        
    elif p['EquationSet']==f_reduced_relax:
        p['Nvar']= 5
        y[0,:] = ic['g']          
        y[1,:] = ic['lambdapoint'] 
        y[2,:] = ic['lambd']       
        y[3,:] = ic['omega']      
        y[4,:] = ic['d']           
            
    elif p['EquationSet']==f_GEMMES:
        p['Nvar'] = 17
        y = np.zeros((p['Nvar'],p['Nx']))
        
        y[ 0,:]= ic['K']               # Capital 
        y[ 1,:]= ic['a']               # Productivity (written A elsewhere)
        y[ 2,:]= ic['N']               # World population
        y[ 3,:]= ic['W']               # Salary 
        y[ 4,:]= ic['P']               # Global price
        y[ 5,:]= ic['D']               # Debt
        
        y[ 6,:]= ic['sigma']           # Intensity of emission in production
        y[ 7,:]= ic['Pbs']             # Price Backstop
        y[ 8,:]= ic['Pc']              # Carbon price
        
        y[ 9,:]= ic['Fexo']            # Exogenous forcing intensity
        y[10,:]= ic['Eland']           # Emission from land    
        y[11,:]= ic['CO2at']           # CO2 in atmosphere
        y[12,:]= ic['CO2up']           # CO2 in upper ocean and biosphere
        y[13,:]= ic['CO2lo']           # CO2 in deep ocean
        y[14,:]= ic['T']               # Atmospheric temperature 
        y[15,:]= ic['T0']              # Deeper ocean temperature 
        y[16,:]= ic['t']               # Time 
        
    elif p['EquationSet']==f_GEMMES_CLIM:
        p['Nvar'] = 8
        y = np.zeros((p['Nvar'],p['Nx']))
      
        y[0,:]= ic['Fexo']            # Exogenous forcing intensity
        y[1,:]= ic['Eland']           # Emission from land    
        y[2,:]= ic['CO2at']           # CO2 in atmosphere
        y[3,:]= ic['CO2up']           # CO2 in upper ocean and biosphere
        y[4,:]= ic['CO2lo']           # CO2 in deep ocean
        y[5,:]= ic['T']               # Atmospheric temperature 
        y[6,:]= ic['T0']              # Deeper ocean temperature 
        y[7,:]= ic['t']               # Time 
                
    elif p['EquationSet']==f_GEMMES_CORE:
        p['Nvar'] = 10
        y = np.zeros((p['Nvar'],p['Nx']))
        
        y[ 0,:]= ic['K']               # Capital 
        y[ 1,:]= ic['a']               # Productivity (written A elsewhere)
        y[ 2,:]= ic['N']               # World population
        y[ 3,:]= ic['W']               # Salary 
        y[ 4,:]= ic['P']               # Global price
        y[ 5,:]= ic['D']               # Debt
        
        y[ 6,:]= ic['sigma']           # Intensity of emission in production
        y[ 7,:]= ic['Pbs']             # Price Backstop
        y[ 8,:]= ic['Pc']              # Carbon price
        
        y[ 9,:]= ic['t']               # Time        
    else : raise Exception("Sorry, No description of this equation set in FG.prepareY !") 
    return y,p
    
###############################################################################
### Get Results ###############################################################
###############################################################################   
def expandY(Y_s,t_s,op,p):
    '''
    TRANFORM THE RAW DATA OF AN EXPERIMENT INTO SOMETHING PRACTICAL.
    LINKS DYNAMIC VARIABLES TO MACRO. 
    
    INPUT : 
        Y_S : Dynamic elements of the simulations
        t_s : time vector of the simulation
        op  : spatial operators
        p   : parameters to control the system. p['EquationSet'] IS PRIMORDIAL
            
    OUTPUT : 
        r   : A dictionnary of all relevant variables 
    '''
    
    '''
    List of variables that can be interesting to have with it : 
        r['Y']      =
        r['K']      =
        r['L']      =
        r['N']      =
        r['g']      =
        r['a']      =
        r['I']      =
        r['W']      =
        r['pi']     =   
        r['Pi']     =
        
        r['omega']  =       
        r['lambda'] =

        #### DEBT
        r['D'] =       # Total debt
        r['d'] =       # Relative debt
        
        #### CES 
        r['nu']=
        r['b'] =
        
        #### PURCHASE POWER
        r['p'] =
        r['i'] =
        
        #### VARIABLES WITH RELAXATION EQUATION
        r['Keff']      =
        r['Leff']      =
        r['lambdaeff'] =
        
        #### LES VARIABLES TYPIQUES DE GEMMES
        r['sigma'] =  # Intensity of emission in production
        r['Pbs']   =  # Price Backstop
        r['Pc']    =  # Carbon price
        r['Fexo']  =  # Exogenous forcing intensity
        r['Eland'] =  # Emission from land    
        r['CO2at'] =  # CO2 in atmosphere
        r['CO2up'] =  # CO2 in upper ocean and biosphere
        r['CO2lo'] =  # CO2 in deep ocean
        r['T']     =  # Atmospheric temperature 
        r['T0']    =  # Deeper ocean temperature     
    ''' 
    
    r= {}
    r['t']=t_s
    
    if p['EquationSet']==f_reduced:
        r['omega'] = Y_s[0,:,:]
        r['lambda']= Y_s[1,:,:]
        r['d']     = Y_s[2,:,:]

        r['nu']    = p['nu']
        r['b'] = .5

        r['i']      = p['eta']*(r['omega']*p['mu']-1)
        r['g']     = (1-r['omega'])/p['nu'] - p['delta1']
        r['pi']    = 1 - r['omega'] - p['r']*r['d']
        
        r['p']      = M.sumexp(r['i']    ,p,r)
        r['Y']      = M.sumexp(r['g']    ,p,r)
        r['a']      = M.sumexp(p['alpha'],p,r)
        r['N']      = M.sumexp(p['beta'] ,p,r)
        
        r['L']      = r['N']*r['lambda']
        r['K']      = r['Y']*r['nu']
        r['W']      = r['a']*r['omega']
        r['Pi']     = r['Y']*r['pi']
        r['D']      = r['Y']*r['d']
        r['I']      = M.fk(r['pi'],p)*r['Y']


        
    elif p['EquationSet']==f_ces_reduced:
        r['omega'] = Y_s[0,:,:]
        r['lambda']= Y_s[1,:,:]
        r['d']     = Y_s[2,:,:]
   
        wp = r['omega']* ( (p['eta']/(1+p['eta']))*(M.philips(r['lambda'],p)-p['alpha']      ))      
        r['g'] = (1-r['omega'])**(1+1/p['nu']) *(1/p['nu'])/(p['b']**(1/p['nu'])) - p['delta1'] - wp/((1-r['omega'])*p['eta'])
        r['pi']    = 1 - r['omega'] - p['r']*r['d']
        r['nu']    = p['nu']* ( (1-r['omega'])/p['b'])**(-1/p['eta'])
        
        r['omega'] = Y_s[0,:,:]
        r['lambda']= Y_s[1,:,:]
        r['d']     = Y_s[2,:,:]
        r['g']     = (1-r['omega'])/p['nu'] - p['delta1']
        r['nu']    = p['nu']* ( (1-r['omega'])/p['b'])**(-1/p['eta'])
        r['pi']    = 1 - r['omega'] - p['r']*r['d']
        
    elif p['EquationSet']==f_reduced_relax:        
        r['g']           = Y_s[0,:,:]
        r['lambdapoint'] = Y_s[1,:,:]
        r['lambd']       = Y_s[2,:,:]
        r['omega']       = Y_s[3,:,:]
        r['d']           = Y_s[4,:,:] 
        
        r['nu']    = p['nu']
        r['pi']    = 1 - r['omega'] - p['r']*r['d']
        
    elif p['EquationSet']==f_full_leontiev:       
        r['a']     = Y_s[0,:,:]
        r['K']     = Y_s[1,:,:]
        r['L']     = Y_s[2,:,:]
        r['W']     = Y_s[3,:,:]
        r['N']     = Y_s[4,:,:]
        r['D']     = Y_s[5,:,:]
        
        r['Yc']    = (1/p['nu'])*r['K']*p['b']**(-1/p['eta'])           # Characteristic value of capital
        r['Lc']    = r['K']/r['a'] * ((1-p['b'])/p['b'])**(1/p['eta'])       # Characteristic quantity of labor
        r['Y']     = r['Yc'] * (1 + (r['L']/r['Lc'])**(-p['eta']))**(-1/p['eta'])  # GDP
        r['omega'] = (r['W']*r['L'])/r['Y']
        r['lambda']= r['L']/r['N']
        r['d']     = r['D']/r['Y']
        #r['pi']    = (r['Y']-r['W']*r['L']-r['D']*p['r'])/r['Y']   
        #r['g']     = np.gradient(r['Y'],p['dt'],axis=0)/r['Y']
        #r['nu']    = r['K']/r['Y']
        
    elif p['EquationSet']==f_full_CES:     
        r['a']     = Y_s[0,:,:]
        r['K']     = Y_s[1,:,:]
        r['L']     = Y_s[2,:,:]
        r['W']     = Y_s[3,:,:]
        r['N']     = Y_s[4,:,:]
        r['D']     = Y_s[5,:,:]
        
        r['Yc'] = r['K']/p['nu']
        r['Lc'] = r['K']/(p['nu']*r['a'])
        r['Y']     = r['Yc'] * (1 + (r['L']/r['Lc'])**(-p['eta']))**(-1/p['eta'])
        r['d']     = r['D']/r['Y']

  
        r['omega'] = (r['W']*r['L'])/r['Y']
        r['lambda']= r['L']/r['N']
        r['d']     = r['D']/r['Y']

        r['pi']    = (r['Y']-r['W']*r['L']-r['D']*p['r'])/r['Y']   
        r['g']     = np.gradient(r['Y'],p['dt'],axis=0)/r['Y']
        r['nu']    = r['K']/r['Y']      
        
    elif p['EquationSet']==f_GEMMES:
        r['K']     = Y_s[ 0,:]          # Capital 
        r['a']     = Y_s[ 1,:]          # Productivity (written A elsewhere)
        r['N']     = Y_s[ 2,:]          # World population
        r['W']     = Y_s[ 3,:]          # Salary 
        r['P']     = Y_s[ 4,:]          # Global price
        r['D']     = Y_s[ 5,:]          # Debt
        
        r['sigma'] = Y_s[ 6,:]          # Intensity of emission in production
        r['Pbs']   = Y_s[ 7,:]          # Price Backstop
        r['Pc']    = Y_s[ 8,:]          # Carbon price
        
        r['Fexo']  = Y_s[ 9,:]          # Exogenous forcing intensity
        r['Eland'] = Y_s[10,:]          # Emission from land    
        r['CO2at'] = Y_s[11,:]          # CO2 in atmosphere
        r['CO2up'] = Y_s[12,:]          # CO2 in upper ocean and biosphere
        r['CO2lo'] = Y_s[13,:]          # CO2 in deep ocean
        r['T']     = Y_s[14,:]          # Atmospheric temperature 
        r['T0']    = Y_s[15,:]          # Deeper ocean temperature 
        r['t']     = Y_s[16,:]          # Time     
        
        r['Y0']    = r['K']/p['nu']
        r['omega'] = (r['W']*r['L'])/r['Y']
        r['lambda']= r['L']/r['N']       
        r['d']     = r['D']/r['Y0']     # CAREFUL HERE WE ARE NOT ON Y ADJUSTED 
        r['g']     = np.gradient(r['Y0'],p['dt'],axis=0)/r['Y0']

    elif p['EquationSet']==f_GEMMES_CLIM:        
        r['Fexo']  = Y_s[0,:]          # Exogenous forcing intensity
        r['Eland'] = Y_s[1,:]          # Emission from land    
        r['CO2at'] = Y_s[2,:]          # CO2 in atmosphere
        r['CO2up'] = Y_s[3,:]          # CO2 in upper ocean and biosphere
        r['CO2lo'] = Y_s[4,:]          # CO2 in deep ocean
        r['T']     = Y_s[5,:]          # Atmospheric temperature 
        r['T0']    = Y_s[6,:]          # Deeper ocean temperature 
        r['t']     = Y_s[7,:]          # Time     
        
    elif p['EquationSet']==f_GEMMES_CORE:
        r['K']     = Y_s[ 0,:]          # Capital 
        r['a']     = Y_s[ 1,:]          # Productivity (written A elsewhere)
        r['N']     = Y_s[ 2,:]          # World population
        r['W']     = Y_s[ 3,:]          # Salary 
        r['P']     = Y_s[ 4,:]          # Global price
        r['D']     = Y_s[ 5,:]          # Debt
        
        r['sigma'] = Y_s[ 6,:]          # Intensity of emission in production
        r['Pbs']   = Y_s[ 7,:]          # Price Backstop
        r['Pc']    = Y_s[ 8,:]          # Carbon price
        
        r['Y0']    = r['K']/p['nu']
        r['omega'] = (r['W']*r['L'])/r['Y']
        r['lambda']= r['L']/r['N']       
        r['d']     = r['D']/r['Y0']     # CAREFUL HERE WE ARE NOT ON Y ADJUSTED 
        r['g']     = np.gradient(r['Y0'],p['dt'],axis=0)/r['Y0']
       
    else : raise Exception("Sorry, wrong Equation set !") 
    return r
