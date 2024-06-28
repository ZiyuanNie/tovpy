#!/usr/bin/env python3

""" 
Various classes for barotropic equations of state (EOS)
"""

import sys, os, shutil
import numpy as np
import scipy as sp
from scipy.interpolate import CubicSpline
import pandas as pd
from io import StringIO
import units

# Data for piecewise polytropes with 4 pieces
PIECEWISE_POLYTROPE_TABLE4 = """
eos logP1 gamma1 gamma2 gamma3
PAL6 34.380 2.227 2.189 2.159 
SLy 34.384 3.005 2.988 2.851 
AP1 33.943 2.442 3.256 2.908 
AP2 34.126 2.643 3.014 2.945
AP3 34.392 3.166 3.573 3.281  
AP4 34.269 2.830 3.445 3.348 
FPS 34.283 2.985 2.863 2.600 
WFF1 34.031 2.519 3.791 3.660 
WFF2 34.233 2.888 3.475 3.517  
WFF3 34.283 3.329 2.952 2.589  
BBB2 34.331 3.418 2.835 2.832 
BPAL12 34.358 2.209 2.201 2.176 
ENG 34.437 3.514 3.130 3.168 
MPA1 34.495 3.446 3.572 2.887 
MS1 34.858 3.224 3.033 1.325 
MS2 34.605 2.447 2.184 1.855 
MS1b 34.855 3.456 3.011 1.425 
PS 34.671 2.216 1.640 2.365 
GS1 34.504 2.350 1.267 2.421 
GS2 34.642 2.519 1.571 2.314 
BGN1H1 34.623 3.258 1.472 2.464 
GNH3 34.648 2.664 2.194 2.304 
H1 34.564 2.595 1.845 1.897
H2 34.617 2.775 1.855 1.858
H3 34.646 2.787 1.951 1.901
H4 34.669 2.909 2.246 2.144
H5 34.609 2.793 1.974 1.915
H6 34.593 2.637 2.121 2.064
H7 34.559 2.621 2.048 2.006
PCL2 34.507 2.554 1.880 1.977 
ALF1 34.055 2.013 3.389 2.033 
ALF2 34.616 4.070 2.411 1.890 
ALF3 34.283 2.883 2.653 1.952 
ALF4 34.314 3.009 3.438 1.803
"""



class EOSPolytropic(object):

    """

    Class for Polytropic EOS

    p = K rho^Gamma

    """
    
    def __init__(self, 
                 polytropic_constant = 1.,
                 polytropic_exponent = 2.):

        self.Gamma = polytropic_exponent
        self.Gammamo = polytropic_exponent - 1.
        self.K = polytropic_constant

    def __give_eos_params(self):
        return self.K, self.Gamma, self.Gammamo

    """ rho-functions """
    
    def Pressure_Of_RestMassDensity(self,rho):
        K,G,G_1 = self.__give_eos_params()
        p = K * rho**G
        return p

    def SpecificInternalEnergy_Of_RestMassDensity(self,rho):
        K,G,G_1 = self.__give_eos_params()
        eps = K/G_1 * rho**G_1 
        return eps

    def EnergyDensity_Of_RestMassDensity(self,rho):
        eps = self.SpecificInternalEnergy_Of_RestMassDensity(rho)
        e = rho * ( 1.0 + eps )
        return e

    def PseudoEnthalpy_Of_BaryonRestMassDensity(self,rho):
        K,G,G_1 = self.__give_eos_params()
        H1 = G*K/G_1 * rho**G_1
        h = np.log(H1+1.0)
        return h

    """ h-functions """
    
    def BaryonRestMassDensity_Of_PseudoEnthalpy(self,h):
        K,G,G_1 = self.__give_eos_params()
        H1 = 2. * np.sinh(0.5*h) * np.exp(0.5*h) # = exp(h) - 1
        rho = ( G_1/(G*K) * H1 )**(1.0/G_1)
        return rho
    
    def Pressure_Of_PseudoEnthalpy(self,h):
        rho = self.BaryonRestMass_Of_PseudoEnthalpy(h)
        p = self.Pressure_Of_RestMassDensity(rho)
        return p

    def SpecificInternalEnergy_Of_PseudoEnthalpy(self,h):
        rho = self.BaryonRestMass_Of_PseudoEnthalpy(h)
        eps = self.SpecificInternalEnergy_Of_RestMassDensity(rho)
        return eps
    
    def EnergyDensity_Of_PseudoEnthalpy(self,h):
        rho = self.BaryonRestMass_Of_PseudoEnthalpy(h)
        e = self.EnergyDensity_Of_RestMassDensity(rho)
        return e

    """ p-functions """
    
    def RestMassDensity_Of_Pressure(self,p):
        K,G,G_1 = self.__give_eos_params()
        rho = (p/K)**(1.0/G)
        return rho
        
    def EnergyDensity_Of_Pressure(self,p):
        rho = self.RestMassDensity_Of_Pressure(p)
        e = self.EnergyDensity_Of_RestMassDensity(rho)
        return e
    
    def EnergyDensityDeriv_Of_Pressure(self,p):
        #TODO
        K,G,G_1 = self.__give_eos_params()
        dedp = (p/K)**(1.0/G)*(1.0+K*(p/K)**(G_1/G)*(p+G_1)/(G_1*p))/G
        return dedp
        
    def PseudoEnthalpy_Of_Pressure(self,p):
        rho = self.RestMassDensity_Of_Pressure(p)
        h = self.PseudoEnthalpy_Of_BaryonRestMassDensity(rho)
        return h 

    # EOSPolytropic ---



class EOSPiecewisePolytropic(object):
    
    """

    Class for Piecewise Polytropic EOS

    """
    
    def __init__(self, name, **params):
        
        self.nPoly = 0 # number of pieces
        self.rhoTab = [] # starting rest-mass density of polytropic piece i (kg/m^3) 
        self.epsTab = [] # starting energy density of polytropic piece i (J/m^3) 
        self.pTab = [] # starting pressure of polytropic piece i (Pa, J/m^3) 
        self.kTab = [] # polytropic constant of piece i 
        self.gammaTab = [] # adiabatic index of piece i 
        self.nTab = [] # polytropic index n_i = 1/(Gamma_i - 1) 
        self.aTab = [] # integration constant (Eq. 3) of PRD 79, 124032 (2009) 
        self.hTab = [] # pseudo-enthalpy

        self.uts = units.Units()
        
        if name == 'piecewise_poly_1' or name == 'poly':

            if "gamma" not in params:
                raise ValueError("The polytropic index 'gamma' is required ")
            if "K" not in params:
                raise ValueError("The polytropic constant 'K' is required")
            
            self.__setup_piecewise_polytrope_1(params['gamma'],params['K'])

        elif name == 'piecewise_poly_4':
            
            if "gamma1" not in params:
                raise ValueError("The polytropic index 'gamma1' is required ")
            if "gamma2" not in params:
                raise ValueError("The polytropic index 'gamma2' is required ")
            if "gamma3" not in params:
                raise ValueError("The polytropic index 'gamma3' is required ")
            if "logP1" not in params:
                raise ValueError("The polytropic constant 'logP1' is required")

            self.__setup_piecewise_polytrope_4(params['logP1'],params['gamma1'],params['gamma2'],params['gamma3'])
            
        else:
                
            params = self.__find_eos_in_existing_piecewise_polytrope_4(name)
            if not params:
                raise ValueError("Unknown EOS name {}".format(name))

            self.__setup_piecewise_polytrope_4(params['logP1'],params['gamma1'],params['gamma2'],params['gamma3'])
        return
    
    """ setup functions """
    
    def __setup_piecewise_polytrope_1(self,gamma,K):
        """
        Single piece polytropic EOS
        """
        self.nPol = 1
        zero = np.atleast_1d(0.)
        one = np.atleast_1d(1.)
        self.rhoTab = zero
        self.epsTab = zero
        self.pTab = zero
        self.kTab = K * one
        self.gammaTab = gamma * one
        self.nTab = one / (gamma - 1.0)
        self.aTab = zero
        self.hTab = zero
        return
    
    def __setup_piecewise_polytrope_4(self,logp1_SI,gamma1,gamma2,gamma3):
        """
        Four pieces polytropic EOS 
        """

        if gamma1 <= 1.0 or gamma2 <= 1.0 or gamma3 <= 1.0:
            raise ValueError("gamma_i > 1")
        
        # Data for the 4-piece piecewise polytrope fitting the low-density part of the SLY4 EOS.  
        # Pressure is defined in Pa (N/m^2) 
        # Rest-mass density is in kg/m^3
        rhoLow = np.array([0, 2.44033979e10, 3.78358138e14, 2.62780487e15])
        kLow = np.array([1.0801158752700761e7, 1.311359898998385e10, 6.507604807550857e19, 3.053461077133694e8])
        gammaLow = np.array([1.58424999, 1.28732904, 0.62223344, 1.35692395])
 
        # Transition densities between the 3 high-density polytropes */
        rho1 = 10.0**17.7
        rho2 = 10.0**18.0
        
        # Pressure at rho1 
        p1 = 10.0**logp1_SI
        
        # Polytropic constants 
        k1 = p1 / rho1**gamma1
        k2 = p1 / rho1**gamma2
        k3 = k2 * rho2**(gamma2 - gamma3)
 
        # Calculate the variable joining density rho0 between the high and low density EOS 
        rho0 = (kLow[3] / k1)**( 1.0 / (gamma1 - gammaLow[3]))
        p1min = kLow[3] * rho1**gammaLow[3]
        
        if logp1_SI < np.log10(p1min) or logp1_SI > 34.5:
            raise ValueError("logp1_SI = {:.3f} should be between {:.3f} and 34.5".format(logp1_SI,np.log10(p1min)))
 
        # Add another polytrope if the joining density is
        # below the start of the last low density polytrope or
        # above the end of the first high density polytrope.
        
        if rho0 > rhoLow[3] and rho0 < rho1:
            
            # No issue. There will be a total of 7 polytropes.
            
            self.nPoly = 7
            #zero = np.zeros(self.nPoly)
            self.kTab = np.zeros(self.nPoly)
            self.gammaTab = np.zeros(self.nPoly)
            self.rhoTab = np.zeros(self.nPoly)
            
            self.kTab[0] = kLow[0]
            self.kTab[1] = kLow[1]
            self.kTab[2] = kLow[2]
            self.kTab[3] = kLow[3]
            self.kTab[4] = k1
            self.kTab[5] = k2
            self.kTab[6] = k3
            
            self.gammaTab[0] = gammaLow[0]
            self.gammaTab[1] = gammaLow[1]
            self.gammaTab[2] = gammaLow[2]
            self.gammaTab[3] = gammaLow[3]
            self.gammaTab[4] = gamma1
            self.gammaTab[5] = gamma2
            self.gammaTab[6] = gamma3
            
            self.rhoTab[0] = rhoLow[0]
            self.rhoTab[1] = rhoLow[1]
            self.rhoTab[2] = rhoLow[2]
            self.rhoTab[3] = rhoLow[3]
            self.rhoTab[4] = rho0
            self.rhoTab[5] = rho1
            self.rhoTab[6] = rho2
            
        else:
            
            # Add an 8th polytrope between gammaLow[3] and gamma1.
            # It will be between the densities rhoJoin1 and rhoJoin2. 

            rhoJoin1 = 5.0e15
            rhoJoin2 = 1.0e16
 
            # Calculate the pressure at the start and end densities. 
            pJoin1 = kLow[3] * rhoJoin1**gammaLow[3]
            pJoin2 = k1 * rhoJoin2**gamma1
 
            # Calculate K and Gamma for the joining polytrope 
            gammaJoin = np.log(pJoin2 / pJoin1) / np.log(rhoJoin2 / rhoJoin1)
            kJoin = pJoin1 / rhoJoin1**gammaJoin

            # Now join all 8 polytropes. 
            self.nPoly = 8
            #zero = np.zeros(self.nPoly)
            self.kTab = np.zeros(self.nPoly) 
            self.gammaTab = np.zeros(self.nPoly)
            self.rhoTab = np.zeros(self.nPoly)

            self.kTab[0] = kLow[0]
            self.kTab[1] = kLow[1]
            self.kTab[2] = kLow[2]
            self.kTab[3] = kLow[3]
            self.kTab[4] = kJoin
            self.kTab[5] = k1
            self.kTab[6] = k2
            self.kTab[7] = k3
            
            self.gammaTab[0] = gammaLow[0]
            self.gammaTab[1] = gammaLow[1]
            self.gammaTab[2] = gammaLow[2]
            self.gammaTab[3] = gammaLow[3]
            self.gammaTab[4] = gammaJoin
            self.gammaTab[5] = gamma1
            self.gammaTab[6] = gamma2
            self.gammaTab[7] = gamma3
            
            self.rhoTab[0] = rhoLow[0]
            self.rhoTab[1] = rhoLow[1]
            self.rhoTab[2] = rhoLow[2]
            self.rhoTab[3] = rhoLow[3]
            self.rhoTab[4] = rhoJoin1
            self.rhoTab[5] = rhoJoin2
            self.rhoTab[6] = rho1
            self.rhoTab[7] = rho2
 
        # Convert to geometric units
        G_SI = self.uts.const('G_SI')
        C_SI = self.uts.const('C_SI')
        G_C2_SI = self.uts.const('C_SI')

        self.rhoTab *= G_C2_SI 
        self.kTab *= G_SI**(1.0 - self.gammaTab) * C_SI**(2.0 * self.gammaTab - 4.0)
        ##print(self.gammaTab)#,self.kTab,self.rhoTab) 
        
        # Calculate remaining quantities (p, n, a, eps, h)
        rho_i = self.rhoTab
        p_i = self.kTab * np.power(rho_i, self.gammaTab)
        n_i = 1.0 / (self.gammaTab - 1.0);
        
        a_i = np.zeros_like(p_i)
        i = np.arange(1,self.nPoly)
        a_i[i] = a_i[i-1] + (n_i[i-1] - n_i[i]) * p_i[i] / rho_i[i]
        
        eps_i = (1.0 + a_i) * rho_i + n_i * p_i
        
        enthalpy_i = 1.0 + a_i + (n_i + 1) * p_i / rho_i
        enthalpy_i[0] = 1.0 # p/rho -> 0 as rho -> 0, and a_0 = 0
        
        self.pTab = p_i
        self.nTab = n_i
        self.aTab = a_i
        self.epsTab = eps_i
        self.hTab = np.log(enthalpy_i)
        
        return
    
    def __existing_piecewise_polytrope_4(self):
        """
        Load data for existing 4-pieces piecewise polytropic
        """
        return pd.read_csv(StringIO(PIECEWISE_POLYTROPE_TABLE4), delim_whitespace=True).to_dict('records')
    
    def __find_eos_in_existing_piecewise_polytrope_4(self, name):
        """
        Load data for existing 4-pieces piecewise polytropic
        """
        pp4data = self.__existing_piecewise_polytrope_4()
        return next((d for d in pp4data if d['eos'] == name), None)
    
    """ h-functions """
    
    def __polytrope_piece_of_h(self,h):
        """
        Determine which polytrope piece h belongs to
        hTab[i] is starting pseudo-enthalpy of polytropic piece i 
        """
        return np.where(h>self.hTab)[0][0]
    
    def RestMassDensity_Of_PseudoEnthalpy(self,h):
        i = self.__polytrope_piece_of_h(h)
        enthalpy = np.exp(h)
        k_i = self.kTab[i]
        n_i = self.nTab[i]
        a_i = self.aTab[i]
        num = enthalpy - 1.0 - a_i
        den = (n_i + 1.0) * k_i
        rho = (num/den)**n_i
        return rho

    def Pressure_Of_PseudoEnthalpy(self,h):
        i = self.__polytrope_piece_of_h(h)
        enthalpy = np.exp(h)
        k_i = self.kTab[i]
        n_i = self.nTab[i]
        a_i = self.aTab[i]
        num = enthalpy - 1.0 - a_i
        den = (n_i + 1.0) * k_i
        p = k_i * (num/den)**(n_i + 1.0)
        return p

    def EnergyDensity_Of_PseudoEnthalpy(self,h):
        i = self.__polytrope_piece_of_h(h)
        enthalpy = np.exp(h)
        k_i = self.kTab[i]
        n_i = self.nTab[i]
        a_i = self.aTab[i]
        rho = self.RestMassDensity_Of_PseudoEnthalpy(h)
        num = 1.0 + a_i + n_i * enthalpy
        den = n_i + 1.0
        e = rho*(num/den)
        return e

    def EnergyDensityDeriv_Of_PseudoEnthalpy(self,h):
        i = self.__polytrope_piece_of_h(h)
        enthalpy = np.exp(h)
        n_i = self.nTab[i]
        a_i = self.aTab[i]
        dedp = (n_i*enthalpy)/(enthalpy - 1.0 - a_i)
        return dedp

    def SoundSpeed_Of_PseudoEnthalpy(self,h):
        i = self.__polytrope_piece_of_h(h)
        enthalpy = np.exp(h)
        n_i = self.nTab[i]
        a_i = self.aTab[i]
        cs = np.sqrt((enthalpy - 1.0 - a_i) / (n_i * enthalpy))
        return cs
    
    """ p-functions """
    
    def __polytrope_piece_of_p(self,p):
        """
        Determine which polytrope piece p belongs to
        pTab[i] is starting pressure of polytropic piece i 
        """
        return np.where(p>self.pTab)[0][0]

    def PseudoEnthalpy_Of_Pressure(self,p):
        i = self.__polytrope_piece_of_p(p)
        k_i = self.kTab[i]
        n_i = self.nTab[i]
        a_i = self.aTab[i]
        enthalpy = 1.0 + a_i + (n_i + 1.0) * k_i * (p / k_i)**(1.0 / (n_i + 1.0))
        h = np.log(enthalpy)
        return h

    def RestMassDensity_Of_Pressure(self,p):
        i = self.__polytrope_piece_of_p(p)
        k_i = self.kTab[i]
        n_i = self.nTab[i]
        rho = (p/k_i)**(n_i/(n_i+1.0))
        return rho

    def EnergyDensity_Of_Pressure(self,p):
        i = self.__polytrope_piece_of_p(p)
        k_i = self.kTab[i]
        n_i = self.nTab[i]
        a_i = self.aTab[i]
        rho = (p/k_i)**(n_i/(n_i+1.0))
        e = (1.0 + a_i) * rho + n_i * p
        return e

    def EnergyDensityDeriv_Of_Pressure(self,p):
        i = self.__polytrope_piece_of_h(h)
        gamma_i = self.gammaTab[i]
        e = self.EnergyDensity_Of_Pressure(p)
        dedp = (epsilon + p) / (gamma_i * p) # returns nan at p=0 
        return dedp

    def SoundSpeed_Of_Pressure(self,p):
        i = self.__polytrope_piece_of_h(h)
        gamma_i = self.gammaTab[i]
        e = self.EnergyDensity_Of_Pressure(p)
        cs = np.sqrt( (gamma_i*p)/(e + p) ) # returns nan at p=0 
        return cs

    """ rho-functions """
    
    def __polytrope_piece_of_rho(self,rho):
        """
        Determine which polytrope piece rho belongs to
        rhoTab[i] is starting rest-mass density of polytropic piece i 
        """
        return np.where(rho>self.rhoTab)[0][0]

    # EOSPiecewisePolytropic ---



class EOSTabular(object):

    """

    Class for tabulated EOS 

    EOS data is expected in the format:

    L1
    C1 C2 C3 C4
    ... 

    with 

    L1 number of tabulated points
    C1 energy density/c^2    (g/cm^3)
    C2 pressure              (dynes/cm^2)
    C3 enthalpy              (cm^2/s^2)
    C4 baryon number density (cm^{-3})

    EOS tables can be found e.g.

    http://xtreme.as.arizona.edu/NeutronStars/index.php/dense-matter-eos/
    https://git.ligo.org/lscsoft/lalsuite/-/tree/master/lalsimulation/lib
    https://github.com/lscsoft/bilby/tree/master/bilby/gw/eos/eos_tables
    https://bitbucket.org/bernuzzi/tov/src/master/EOS/

    #TODO adapt code for this format 
    #     adapt LAL tables to this format

    """

    def __init__(self, name, **params):

        self.table = []
        self.uts = units.Units()

        if name == 'from_file':

            if "filename" not in params:
                raise ValueError("Need a filename")
            
            self.table = np.loadtxt(params['filename'],
                                    skiprows = 1)
            
        elif name == 'from_ndarray':

            if "data" not in params:
                raise ValueError("Need data")

            self.table = params['data']
            
        else:
            
            raise ValueError("Unknown EOS name")
        
        self.table = self.__remove_leading_zero(self.table)

        self.min_pTab = []
        self.min_eTab = []
        self.min_hTab = []

        self.interp_EnergyDensity_from_Pressure = []
        self.interp_EnergyDensity_from_PseudoEnthalpy = []
        self.interp_Pressure_from_PseudoEnthalpy = []
        self.interp_PseudoEnthalpy_from_EnergyDensity = []

        self.__interpolate()
        
    def __interpolate(self,kind='cubic'):
        """
        Uses a CubicSpline, for which derivatives are also provided
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html
        """
        
        eTab = self.table[:, 0]
        pTab = self.table[:, 1]
        if self.table.shape[2] == 3:
            hTab = self.__pseudoenthalpy_from_p_and_e(pTab,eTab)
        # ...
        
        self.min_pTab = min(pTab)
        self.min_eTab = min(eTab)
        self.min_hTab = min(hTab)

        self.logpTab = np.log(pTab)
        self.logeTab = np.log(eTab)
        self.loghTab = np.log(hTab)

        self.interp_EnergyDensity_from_Pressure = CubicSpline(logpTab,logeTab)
        self.interp_EnergyDensity_from_PseudoEnthalpy = CubicSpline(loghTab,logeTab)
        self.interp_Pressure_from_PseudoEnthalpy = CubicSpline(loghTab,logpTab)
        self.interp_PseudoEnthalpy_from_EnergyDensity = CubicSpline(logeTab,loghTab)

    def __remove_leading_zero(self, table):
        """
        For interpolation of LALSimulation tables;
        loglog interpolation breaks if the first entries are 0s
        """
        if table[0, 0] == 0. or table[0, 1] == 0.:
            return table[1:, :]
        else:
            return table

    def __pseudoenthalpy_from_p_and_e(self,p,e):
        integrand = p / (e + p)
        return np.cumtrapz(integrand, np.log(p), initial=0) + integrand[0]
        
    """ p-functions """
        
    def EnergyDensity_Of_Pressure(self,p):
        """
        Use non-relativistic degenerate gas, p = K * e**(5./3.)
        and return e = K * p**(3./5.) below min pressure        
        """
        p = np.atleast_1d(p)
        e = np.zeros(p.size)
        ids_lt_min = np.nonzero(np.logical_and(p < self.min_pTab, p > 0.0))
        ids_gt_min = np.nonzero(p >= self.min_pTab)
        e[ids_lr_min] = np.exp( self.logeTab[0] + (3. / 5.) *(self.logpTab[ids_lt_min] - self.logpTab[0]) )
        e[ids_gt_min] = np.exp( self.interp_EnergyDensity_from_Pressure(np.log(p[ids_gt_min])) )
        if e.size == 1:
            return e[0]
        return e
                        
    def EnergyDensityDeriv_Of_Pressure(self,p,finite_diff=False):
        """
        Use non-relativistic degenerate gas, p = K * e**(5./3.)     
        and return drvt of e = K * p**(3./5.) below min pressure    
        """
        if finite_diff:
            return self.EnergyDensityDeriv_Of_Pressure_fd(p)
        p = np.atleast_1d(p)
        dedp = np.zeros(p.size)
        ids_lt_min = np.nonzero(np.logical_and(p < self.min_pTab, p > 0.0))
        ids_gt_min = np.nonzero(p >= self.min_pTab)
        dedp[ids_lr_min] = (3.0 / 5.0) * np.exp(self.logeTab[0] - self.logpTab[0]);
        logp = np.log(p[ids_gt_min])
        loge = self.interp_EnergyDensity_from_Pressure(logp)
        dloge_dlogp = self.interp_EnergyDensity_from_Pressure(logp, 1) # drvt
        dedp[ids_gt_min] = dloge_dlogp * npexp(loge - logp)
        if dedp.size == 1:
            return dedp[0]
        return dedp

    def EnergyDensityDeriv_Of_Pressure_fd(self,p):
        """
        A finite difference method implemented in e.g. bilby
        """
        rel_dp = 1e-5
        dp = p * rel_dp
        e_upper = self.Energy_Of_Pressure(p + dp)
        e_lower = self.Energy_Of_Pressure(p - dp)
        dedp = (eps_upper - eps_lower) / (2. * dp)
        return dedp
        
    def PseudoEnthalpy_Of_Pressure(self,p):
        """
        Use non-relativistic degenerate gas, h = K * p**(2./5.), below min p
        """
        p = np.atleast_1d(p)
        h = np.zeros(p.size)
        ids_lt_min = np.nonzero(np.logical_and(p < self.min_pTab, p > 0.0))
        ids_gt_min = np.nonzero(p >= self.min_pTab)
        logp = log(p)
        h[ids_lt_min] = np.exp(self.loghTab[0] + 0.4 * (logp - self.logpTab[0]))
        logh = self.interp_PseudoEnthalpy_from_Pressure(logp)
        h[ids_gt_min] = np.exp(logh)
        if h.size == 1:
            return h[0]
        return h

    """ h-functions """
    
    def EnergyDensity_Of_PseudoEnthalpy(self,h):
        """
        Use non-relativistic degenerate gas, e = K * h**(3./2.), below min h
        """
        h = np.atleast_1d(h)
        e = np.zeros(h.size)
        ids_lt_min = np.nonzero(np.logical_and(h < self.min_hTab, h > 0.0))
        ids_gt_min = np.nonzero(h >= self.min_hTab)
        logh = np.log(h)
        e[ids_lt_min] = np.exp(self.logeTab[0] + 1.5 * (logh[ids_lt_min] - loghTab[0]))
        loge = self.interp_EnergyDensity_from_PseudoEnthalpy(logh)
        e[ids_gt_min] = np.log(loge)
        if e.size == 1:
            return e[0]
        return e
            
    def Pressure_Of_PseudoEnthalpy(self,h):
        """
        Use non-relativistic degenerate gas, p = K * h**(5./2.), below min h
        """
        h = np.atleast_1d(h)
        p = np.zeros(h.size)
        ids_lt_min = np.nonzero(np.logical_and(h < self.min_hTab, h > 0.0))
        ids_gt_min = np.nonzero(h >= self.min_hTab)
        log_h = log(h)
        p[ids_lt_min] = np.exp(self.logpTab[0] + 2.5 * (logh - self.loghTab[0]))
        logp = self.interp_Pressure_from_PseudoEnthalpy(logh)
        p[ids_gt_min] = exp(logp)
        if p.size == 1:
            return p[0]
        return p
    
    def RestMassDensity_Of_PseudoEnthalpy(self,h):
        """
        Use non-relativistic degenerate gas, rho = K * h**(3./2.), below min h
        #TODO unfinished
        """
        h = np.atleast_1d(h)
        rho = np.zeros(h.size)
        ids_lt_min = np.nonzero(np.logical_and(h < self.min_hTab, h > 0.0))
        ids_gt_min = np.nonzero(h >= self.min_hTab)
        logh = np.log(h)
        #rho[ids_lt_min] = np.exp(self.logrhoTab[0] + 1.5 * (logh - self.loghTab[0]))
        #logrho = #self.interp_RestMassDensity_from_PseudoEnthalpy(logh)
        #rho[ids_gt_min] = np.log(logrho)
        if rho.size == 1:
            return rho[0]
        return rho

    # EOSTabular ---



class EOS(object):
    
    """

    Class wrapping various EOS classes 

    Contains the calls needed in TOVSolver()

    """

    def __init__(self, name, **params):

        self.eos = []
        
        if name in ['poly','piecewise_poly_1','piecewise_poly_4']:
            self.eos = EOSPiecewisePolytropic(name,params)
        elif name in ['tabular']:
            self.eos = EOSTabular(name,params)
        else:
            raise ValueError("Unknown EOS name {}".format(name))

    """ p-functions """
        
    def EnergyDensity_Of_Pressure(self,x):
        return self.eos.EnergyDensity_Of_Pressure(x)
                        
    def EnergyDensityDeriv_Of_Pressure(self,x):
        return self.eos.EnergyDensityDeriv_Of_Pressure(x)

    def PseudoEnthalpy_Of_Pressure(self,x):
        return self.eos.PseudoEnthalpy_Of_Pressure(x)

    """ h-functions """
    
    def EnergyDensity_Of_PseudoEnthalpy(self,x):
        return self.eos.EnergyDensity_Of_PseudoEnthalpy(x)
            
    def Pressure_Of_PseudoEnthalpy(self,x):
        return self.eos.Pressure_Of_PseudoEnthalpy(x)

    # EOS ---



if __name__ == "__main__":


    #TODO extensive tests
    
    eos = EOSPiecewisePolytropic('piecewise_poly_1',gamma=2,K=100)
    #eos = EOSPiecewisePolytropic('piecewise_poly_1',**{'gamma':2,'K':100})
    print(eos.EnergyDensity_Of_Pressure(1e-3))

    ##eos = EOSPiecewisePolytropic('SLy')
    ##print(eos.EnergyDensity_Of_Pressure(1e-3))

