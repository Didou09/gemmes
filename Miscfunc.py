# -*- coding: utf-8 -*-
import numpy as np
from scipy.sparse import diags
import shutil
from datetime import datetime
import os
import pickle
import FunctionsGoodwin as FG 

###############################################################################
### SYSTEM INITIALISATION ###############################################
###############################################################################

def preparePhilips(p):
    """
    Calculate the parameters of philips and solow point
    """
    if p['phiType'] == 'linear' : 
        p['phi0']       = p['phinul']*p['phislope']                # Based on Phinul value
        p['phi1']       = p['phislope']                            # Based on Phinul value
        p['lambdamin']  = (p['phi0']+p['alpha'])/p['phi1'] # Same                
    elif p['phiType'] == 'divergent' : 
        p['phi0']       = p['phinul']    / (1- p['phinul']**p['phislope'])     # Based on Phinul value
        p['phi1']       = p['phinul']**(p['phislope']+1) / (1- p['phinul']**p['phislope'])     # Based on Phinul value
        p['lambdamin']  = 1- np.sqrt( p['phi1']/(p['alpha']+p['phi0'])) # Same
    else : raise Exception("Sorry, wrong Philips-type name in FG.preparePhilips") 
    p['omega0']     = 1-p['nu']*(p['alpha']+p['beta']+p['delta1'])  # Solow point of a classic goodwin system        
    return p

def preparedT(p):
    """
    determine best dt value : for each phenomena, it gets the "typical time", and use it a a basis
    Then we divide this typical time by a safety factor (10 is enough)
    And we calculate the number of iterations needed
    """
    p['dtmax']      = 1/np.amax((np.amax(p['alpha']),np.amax(p['beta']),np.amax(p['delta1']),np.amax(philips(p['lambdamax'],p)) ))
    p['dtmax']     /= p['time_safety']
    
    p['Nt']         = int(p['Tmax']/p['dtmax'])                     # Number of temporal iteration
    p['dt']         = p['dtmax']                                    # We have two conventions...
    p['Ns']         = int(p['Tmax']/p['Tstore'])+1                  # Number of elements stored
    if p['StorageMode'] == 'Full' : 
        p['Ns']    = p['Nt']
        p['Tstore']= p['dt']
    return p


def prepareOperators(p):
    """
    Creation of few spatial operators :
    *    3diag is a local meaning on both neighbors and centered value 
    *    diff_cent is a centered difference operator 
    *    laplacian is a classical laplacian operator
    """
    operator={}
    if p['Nx']>=3:
        operator['3diag']    =diags([1/3*np.ones(p['Nx']-1), 1/3*np.ones(p['Nx']-1),np.ones(p['Nx'])*1/3,1/3,1/3], [1, -1, 0,p['Nx'],-p['Nx']])    
        operator['diff_cent']=diags([0.5*np.ones(p['Nx']-1),-0.5*np.ones(p['Nx']-1),-0.5,0.5]                 , [1, -1   ,p['Nx'],-p['Nx']])
        operator['laplacian']=diags([1*np.ones(p['Nx']-1), -1*np.ones(p['Nx']-1),-2*np.ones(p['Nx'])*1/3,1,1], [1, -1, 0,p['Nx'],-p['Nx']])
        operator['Mean']  = np.ones((p['Nx'],p['Nx']))/(p['Nx'])
    
    else : 
        #operator['3diag']    =diags([1/3*np.ones(p['Nx']-1), 1/3*np.ones(p['Nx']-1),np.ones(p['Nx'])*1/3,1/3,1/3], [1, -1, 0,p['Nx'],-p['Nx']])    
        #operator['diff_cent']=diags([0.5*np.ones(p['Nx']-1),-0.5*np.ones(p['Nx']-1),-0.5,0.5]                 , [1, -1   ,p['Nx'],-p['Nx']])
        #operator['laplacian']=diags([1*np.ones(p['Nx']-1), -1*np.ones(p['Nx']-1),-2*np.ones(p['Nx'])*1/3,1,1], [1, -1, 0,p['Nx'],-p['Nx']])
        operator['Mean']  = np.ones((p['Nx'],p['Nx']))/(p['Nx'])
            

    return operator


###############################################################################
### SEPARATE BEHAVIOR FUNCTIONS ###############################################
###############################################################################

def philips(lamb,p):
    """
    lambda the evolution of the salary, philips(lambda) the rate of salary evolution
    it diverges for lambda=1, lambda \in [0,1[ 
    
    relevant parameters :
        * p['phi0']   : the value when there is no employement
        * p['phi1']   : "sensibility" of employement
    """
    if p['phiType'] == 'linear' :      return - p['phi0'] + p['phi1']
    elif p['phiType'] == 'divergent' : return - p['phi0'] + p['phi1']/ (1-lamb)**(p['phislope']) 

def fk(pi,p):
    """
    Investment function according to the debt 
    When r=0, the system behave as a pure Keen (bypass of kappa). 
    BE CAREFUL IF R endogeneisation !
    """
    
    if p['r']!=0 : return p['k0']+p['k1']*np.exp(p['k2']*pi)
    else :         return pi

def delta(pi,p):
    """
    Portion of profit going into investor
    """
    print('DELTA FUNCTION NOT CODED')
    return p['div0']+pi*p['div1']

###############################################################################
### NUMERICAL CORE ############################################################
###############################################################################
def rk4(y,op,p):
    """
    a traditional RK4 scheme, with y the vector values, and p the parameter dictionnary
    dt is contained within p
    """
    dy1 =  p['EquationSet'](y       ,op,p  )
    dy2 =  p['EquationSet'](y+dy1/2 ,op,p  )
    dy3 =  p['EquationSet'](y+dy2/2 ,op,p  )
    dy4 =  p['EquationSet'](y+dy3   ,op,p  )
    return (dy1 + 2*dy2 + 2*dy3 + dy4)/6 

def sumexp(f,p,r):
    y=np.zeros((p['Nt'],p['Nx']))
    y[0]=1
   
    
    if type(f)==float :
        for i in range(p['Nt']-1): y[i+1,:]=y[i,:]*(1+f   *(r['t'][i+1]-r['t'][i]))      
    else :
        for i in range(p['Nt']-1): y[i+1,:]=y[i,:]+(1+f[i]*(r['t'][i+1]-r['t'][i]))
    return y
###############################################################################
################## MISCELLANEOUS ##############################################
###############################################################################
def getperiods(r,p,op):
    '''
    calculate period index values (index of the local maximal position on lambda)
    '''
    r['PeriodID']=[]
    
    for j in np.arange(p['Nx']): 
        r['PeriodID'].append([])
        id1=1
        while id1<p['Nt']-2:
            if ( r['lambda'][id1,j]>r['lambda'][id1-1,j] and 
                 r['lambda'][id1,j]>r['lambda'][id1+1,j] ):
                r['PeriodID'][j].append(1*id1)  
            id1+=1    
    return r        
            


def ShowtimeParameters(p):
    print(50*'#')
    print('Typical characteristic times of the system in years :')
    if p['alpha'] :print('Population growth        :', 1/p['alpha'])
    if p['beta'] :print('Labour efficiency        :', 1/p['beta'])
    print('Firing                   :', p['tauF'])
    print('Recruitment              :', p['tauR'])
    if p['delta1'] :print('Capital destruction      :', 1/p['delta1'])
    if p['r']     :print('Debt growth              :', 1/p['r'])
    print('Salary fastest inflation :', 1/philips(p['lambdamax'],p))
    print(30*'#')
    print('Chosen dt                :',p['dtmax'])      
    print('Number of step           :',p['Nt'])
    print('Number of systems        :',p['Nx'])
    print('Simulated time           :',p['Tmax'],'years ')
    print('Storage every            :',p['Tstore'])
    print(50*'#')

def ShowSpecificity(p):
    print(50*'#')
    #print('Productive function :',p['GDPtype'])  
    print('Philips function    :',p['phiType'])
    print('System solved       :',end='')
    
    if p['EquationSet']   == FG.f_reduced : print(' GoodwinKeen-Leontiev-reduced')
    elif p['EquationSet'] == FG.f_ces_reduced : print(' GoodwinKeen-CES-reduced')
    elif p['EquationSet'] == FG.f_full_leontiev : print('Full-Leontiev')
    elif p['EquationSet'] == FG.f_full_CES : print('Full-CES')


    print(50*'#')
    print('Migration sensibility  :', p['muN'])    
    print('Investment sensibility :', p['muI'])    
    print('Screened Philips       :', p['g2'] )    
    print('Global Philips         :', p['g1'] )      
     
    
def savedata(rootfold,t,Y_s,p,op):  
    date = str(datetime.now()).split(' ')[0]
    try : os.mkdir(rootfold+date)
    except BaseException : print('Folder of the day already created')
    Nfolder = os.listdir(rootfold+date+'/')
    Fold    = rootfold+date+'/Expe_'+str(len(Nfolder))
    os.mkdir(Fold)   
    pickle.dump(t  ,      open(Fold+'/t.p'      ,"wb"))
    pickle.dump(Y_s,      open(Fold+'/Y_s.p'    ,"wb"))
    pickle.dump(p  ,      open(Fold+'/params.p' ,"wb"))
    pickle.dump(op ,      open(Fold+'/spatialop.p',"wb"))
    
    codefolder = os.getcwd()
    shutil.copyfile(codefolder+'/FunctionsGoodwin.py'     ,Fold+'/FunctionsGoodwin.py' )
    shutil.copyfile(codefolder+'/Main.py'                 ,Fold+'/Main.py' )
    print('data saved in :',Fold)


