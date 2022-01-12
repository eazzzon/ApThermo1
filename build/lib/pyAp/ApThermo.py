# ApThermo model of Li and Costa (2020, GCA)
import scipy.optimize
import math
import numpy as np

## constants
MassCl = 35.45;   #molar mass of Cl
MassF = 18.998;   #molar mass of Cl
MassH2O = 18.015; #molar mass of h2o
meanM = 33;       #molar mass of studied melt
R = 8.314;        #gas constant

Wg_ClOH=5;        # in kJ/mol
Wg_FOH=7;         # in kJ/mol
Wg_FCl=16;        # in kJ/mol


def conversion(moleOH_melt,k2):
    """

    """

    def func_speci(x):
        y = 2*x+(8*x+k2-2*x*k2-math.sqrt(k2)*math.sqrt(16*x-16*x**2+k2-4*x*k2+4*x**2*k2))/(k2-4) - moleOH_melt
        return y

    try:
        moleH2O_melt = scipy.optimize.fsolve(func_speci,0.1)
    except:
        moleH2O_melt = [np.nan, np.nan]

    def func_wt(x):
        y = (x/MassH2O)/(x/MassH2O+(1-x)/meanM) -  moleH2O_melt
        return y
    try:
        MeltWater = scipy.optimize.fsolve(func_wt,0.1)*100     # x100 wt.%
    except:
        MeltWater = [np.nan, np.nan]
    
    return MeltWater[0]


def MeltWater(x_f,x_cl,MeltF,MeltCl,T,kEq):
    
    # calculate Gibbs free energy of reaction
    T_K = T + 273.15;                         # in Kelvin
    deltaG_ClOH = 72.9 - 0.034 * T_K;           # in kJ/mol
    deltaG_FOH = 94.6 - 0.04 * T_K;            # in kJ/mol
    deltaG_ClF = deltaG_FOH - deltaG_ClOH;     # in kJ/mol
  
    Wg_diff1 = Wg_FOH-Wg_FCl;           
    Wg_diff2 = Wg_ClOH-Wg_FCl;           
    Wg_diff3 = Wg_ClOH-Wg_FOH;          
    

    # calculate mole fraction of OH using Cl and F
    x_oh = 1 - (x_cl + x_f)
    if x_oh < 0:
        x_oh = 0
    
    # calculate Kds
    Kd_OHCl= math.exp((1000*(-deltaG_ClOH-((x_cl-x_oh)*Wg_ClOH+x_f*Wg_diff1)))/(R*T_K));
    Kd_OHF = math.exp(1000*(-deltaG_FOH-((x_f-x_oh)*Wg_FOH+x_cl*Wg_diff2))/(R*T_K));
    Kd_ClF = Kd_OHF/Kd_OHCl;
    

    # calculated activity coefficients (gamma)
    gammaOH=math.exp(1000*((x_cl*(1-x_oh)*Wg_ClOH+x_f*(1-x_oh)*Wg_FOH-x_cl*x_f*Wg_FCl))/(R*T_K))
    gammaF=math.exp(1000*((x_cl*(1-x_f)*Wg_FCl+x_oh*(1-x_f)*Wg_FOH-x_cl*x_oh*Wg_ClOH))/(R*T_K))
    gammaCl=math.exp(1000*((x_oh*(1-x_cl)*Wg_ClOH+x_f*(1-x_cl)*Wg_FCl-x_f*x_oh*Wg_FOH))/(R*T_K))
    

    # calculate molar ratios of volatiles in the melt
    ClF_melt = (x_cl/x_f)/Kd_ClF;          
    massClF    = ClF_melt*MassCl/MassF
 
    OHCl_melt=(x_oh/x_cl)/Kd_OHCl  
    OHF_melt=(x_oh/x_f)/Kd_OHF 


    k2 = kEq
    
    if MeltCl == MeltCl:
        moleCl_melt =((MeltCl/10000)/MassCl)/(100/meanM)  # mole Cl in melt (input Cl in ppm)
        moleOH_melt1 = moleCl_melt*OHCl_melt
        MeltWater_Cl = conversion(moleOH_melt1,k2)
    else:
        MeltWater_Cl = float('nan')
    
    if MeltF == MeltF:
        moleF_melt =((MeltF/10000)/MassF)/(100/meanM);    # mole F in melt (input F in ppm)
        moleOH_melt2 = moleF_melt*OHF_melt;
        MeltWater_F = conversion(moleOH_melt2,k2)
    else:
        MeltWater_F = float('nan')
    
    return MeltWater_F,MeltWater_Cl    # x100 wt.%



def ErrorF(x_f,x_cl,MeltF,T,kEq,p1,p2):
    
    # Gibbs free energy of reaction 
    T_K = T + 273.15;                         # in Kelvin
    
    deltaG_FOH = (94.6 + p1) - (0.040 + p2) *T_K;            # in kJ/mol   
  
    Wg_diff1 = Wg_FOH-Wg_FCl;           
    Wg_diff2 = Wg_ClOH-Wg_FCl;           
    Wg_diff3 = Wg_ClOH-Wg_FOH;          
    
    x_oh = 1 - (x_cl + x_f)
    if x_oh < 0:
        x_oh = 0
    
    Kd_OHF = math.exp(1000*(-deltaG_FOH-((x_f-x_oh)*Wg_FOH+x_cl*Wg_diff2))/(R*T_K));
    
    # calculate molar ratios of volatiles in the melt
    OHF_melt=(x_oh/x_f)/Kd_OHF 
    moleF_melt =((MeltF/10000)/MassF)/(100/meanM);    # mole F in melt (input F in ppm)
    moleOH_melt2 = moleF_melt*OHF_melt;
    
    k2 = kEq

    MeltWater_F = conversion(moleOH_melt2,k2)
    
    return MeltWater_F    # x100 wt.%


def ErrorCl(x_f,x_cl,MeltCl,T,kEq,p1,p2):
    
    # Gibbs free energy of reaction 
    T_K = T + 273.15;                         # in Kelvin
    
    deltaG_ClOH     = (72.9 + p1) - (0.034 + p2) *T_K
 
    Wg_diff1 = Wg_FOH-Wg_FCl;           
    Wg_diff2 = Wg_ClOH-Wg_FCl;           
    Wg_diff3 = Wg_ClOH-Wg_FOH;          
    
    x_oh = 1 - (x_cl + x_f)
    if x_oh < 0:
        x_oh = 0
    
    Kd_OHCl= math.exp((1000*(-deltaG_ClOH-((x_cl-x_oh)*Wg_ClOH+x_f*Wg_diff1)))/(R*T_K))
    
    OHCl_melt=(x_oh/x_cl)/Kd_OHCl  
                                            
    moleCl_melt =((MeltCl/10000)/MassCl)/(100/meanM)  # mole Cl in melt (input Cl in ppm)
    moleOH_melt1 = moleCl_melt*OHCl_melt
    
    k2 = kEq

    MeltWater_Cl = conversion(moleOH_melt1,k2)
    
    return MeltWater_Cl    # x100 wt.%
