#!/usr/bin/env python3

""" 
Various classes for barotropic equations of state (EOS)
"""

import sys, os, shutil
import numpy as np
import scipy as sp
from scipy.integrate import solve_ivp
from scipy.special import factorial2
import matplotlib.pyplot as plt

import eos as EOS
import units




class TOV(object):
    
    """ 
    
    Class to solve the Tolman-Oppenheimer-Volkov stellar structure
    equations together with even/odd parity stationary bartropic perturbations
    
    Lindblom , Astrophys. J. 398 569. (1992) 
    Damour & Nagar, Phys. Rev. D 80, 084035 (2009)
    
    Work in geometric units

    Reference codes:

    https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_neutron_star_t_o_v_8c_source.html
    https://bitbucket.org/bernuzzi/tov/src/master/TOVL.m
    https://lscsoft.docs.ligo.org/bilby/_modules/bilby/gw/eos/tov_solver.html
    
    """

    def __init__(self,
                 eos = None, # EOS instance
                 leven = [], # multipole indexes of even perturbations 
                 lodd = [], # multipole indexes of odd perturbations 
                 dhfact = -1e-3, # ODE step
                 ode_method = 'DOP853',
                 ode_atol = 0.,
                 ode_rtol = 1e-9, 
                 output = []): 

        if not eos:
            raise ValueError("Must provide a EOS")
        self.eos = eos

        # Solve perturbation equatins for these indexes
        self.leven = leven[leven>1]
        self.lodd = lodd[lodd>1]

        # Build variable list
        var = self.__buildvars()
        self.nvar = len(var)
        self.var = dict(zip(v,range(nvar)))
        self.ivar = {v: k for k, v in self.var.items()}

        # ODE solver options
        if dhfact > 0.:
            raise ValueError("ODE timestep must be negative")
        self.dhfact = dhfact
        self.ode_method = ode_method        
        self.ode_atol = ode_atol
        self.ode_rtol = ode_rtol

        # output dir
        self.output = output 
        self.__output_makedir(self,outdir=output)

        
    def __buildvars(self):
        """
        List of varnames
        """
        v = ['r','m','nu']
        for l in self.leven:
            v.append('H{}'.format(l))
            v.append('dH{}'.format(l))
        for l in self.lodd:
            v.append('Psi{}'.format(l))
            v.append('dPsi{}'.format(l))
        return v
    
    def __pert_even(self,ell,m,r,p,e,dedp,dnu_dr=[]):
        """
        Eq.(27-29) of Damour & Nagar, Phys. Rev. D 80, 084035 (2009)
        https://arxiv.org/abs/0906.0096
        Note only C0 depends on ell: return an array of values
        """
        r2 = r**2
        r3 = r * r2
        div_r = 1.0/r
        div_r2 = div_r**2
        exp_lam = 1.0 / (1.0 - 2.0 * m * div_r ) # Eq. (18)
        if not dnu_dr:
            dnu2 = (2.0 * (m + 4.0 * np.pi * r3 * p) / (r * (r - 2.0 * m)))**2
        else:
            dnu2 = dnu_dr**2
        C1 = 2.0/r + exp_lam * ( 2*m*div_r2 + 4*np.pi*r*(p-e) ) 
        C0 = np.zeros_like(np.amax(ell))
        for l in ell:
            Lam = l*(l+1)
            C0[l] = -dnu2
            C0[l] += exp_lam * ( -Lam*div_r2 + 4*np.pi*( 5*e + 9*p + (e + p) * dedp ) ) 
        return C1, C0
                
    def __pert_odd(self,ell,m,r,p,e,dedp):
        """
        Eq.(31) of Damour & Nagar, Phys. Rev. D 80, 084035 (2009)
        https://arxiv.org/abs/0906.0096
        Note only C0 depends on ell: return an array of values
        """
        r2 = r**2
        r3 = r * r2
        div_r = 1.0/r
        div_r2 = div_r**2
        div_r3 = div_r*div_r2
        exp_lam = 1.0 / (1.0 - 2.0 * m * div_r ) # Eq. (18)
        C1 = exp_lam * ( 2*m + 4*np.pi*r3*(p-e) ) * div_r2
        C0 = np.zeros_like(np.amax(ell))
        for l in ell:
            Lam = l*(l+1)
            C0[l] = exp_lam*( -Lam*div_r2 + 6*m*div_r3 - 4*np.pi*(e-p) )
    
    def __tov_rhs(self,x,y):
        """
        ODE r.h.s. for TOV equations with pseudo-enthalpy independent variable.
        Implements Eqs. (5) and (6) of Lindblom, Astrophys. J. 398, 569 (1992).
        Also uses Eqs. (7) and (8) [ibid] for inner boundary data, and
        Eqs. (18), (27), (28) of Damour & Nagar, Phys. Rev. D 80, 084035 (2009)
        for the metric perturbation used to obtain the Love number. 
        """
        dy = np.zeros_like(y)
        # Unpack y
        r = y[self.ivar['r']]
        m = y[self.ivar['m']]
        nu = y[self.ivar['nu']]
        # EOS call
        p = self.eos.PressureOfPseudoEnthalpyGeometerized(h)
        e = self.eos.EnergyDensityOfPseudoEnthalpyGeometerized(h)
        dedp = self.eos.EnergyDensityDerivOfPressureGeometerized(p);
        # TOV
        dr_dh = -r * (r - 2.0 * m) / (m + 4.0 * np.pi * r**3 * p)
        dm_dh = 4.0 * np.pi * r**2 * e * dr
        dnu_dr =  2.0 * (m + 4.0 * np.pi * r**3 * p) / (r * (r - 2.0 * m))
        dy[self.ivar['r']] = dr_dh
        dy[self.ivar['m']] = dm_dh
        dy[self.ivar['nu']] = dnu_dr * dr_dh        
        # Even perturbations
        C0,C1 = self.__pert_even(self.leven,m,r,p,e,dedp,dnu_dr)
        for l in self.leven:
            H = y[self.ivar['H{}'.format(l)]]
            dH = y[self.ivar['dH{}'.format(l)]]
            dH_dh = dH * dr_dh;
            ddH_dh = -(C0[l] * H + C1 * dH) * dr_dh;
            dy[self.ivar['H{}'.format(l)]] = dH_dh
            dy[self.ivar['dH{}'.format(l)]] = ddH_dh    
        # Odd perturbations
        C0,C1 = self.__pert_odd(self.leven,m,r,p,e,dedp)
        for l in self.lodd:
            Psi = y[self.ivar['Psi{}'.format(l)]]
            dPsi = y[self.ivar['dPsi{}'.format(l)]]
            dPsi_dh = dH * dr_dh;
            ddPsi_dh = -(C0[l] * Psi + C1 * dPsi) * dr_dh;
            dy[self.ivar['Psi{}'.format(l)]] = dPsi_dh
            dy[self.ivar['dPsi{}'.format(l)]] = ddPsi_dh    
        return dy

    def __initial_data(self,pc,dhfact=-1e-3,verbose=False):
        """
        Set initial data for the solution of TOV equations using the pseudo-enthalpy formalism introduced in:
        Lindblom (1992) "Determining the Nuclear Equation of State from Neutron-Star Masses and Radii", Astrophys. J. 398 569.
        * input the central pressure
        """
        y = np.zeros(self.nvar)
        # Central values 
        ec = self.eos.EnergyDensityOfPressureGeometerized(pc)
        hc = self.eos.PseudoEnthalpyOfPressureGeometerized(pc)
        dedp_c = self.eos.EnergyDensityDerivOfPressureGeometerized(pc)
        dhdp_c = 1.0 / (ec + pc)
        dedh_c = dedp_c / dhdp_c
        dh = dh_fact * hc
        h0 = hc + dh
        h1 = 0.0 - dh
        r0 = np.sqrt(-3.0 * dh / (2.0 * np.pi * (ec + 3.0 * pc)))
        m0 = 4.0 * np.pi * r0**3 * ec / 3.0
        # Series expansion for the initial core 
        r0 *= 1.0 + 0.25 * dh * (ec - 3.0 * pc  - 0.6 * dedh_c) / (ec + 3.0 * pc) # second factor Eq. (7) of Lindblom (1992) 
        m0 *= 1.0 + 0.6 * dh * dedh_c / ec # second factor of Eq. (8) of Lindblom (1992) 
        y['r'] = r0
        y['m'] = m0
        y['nu'] = 0.0
        #  Initial data for the ell-perturbation
        a0 = 1.0
        for l in self.leven:
            y[self.ivar['H{}'.format(l)]] = a0 * r0**l
            y[self.ivar['dH{}'.format(l)]] = a0 * l * r0**(l-1)
        for l in self.leven:
            y[self.ivar['Psi{}'.format(l)]] = a0 * r0**(l+1)
            y[self.ivar['dPsi{}'.format(l)]] = a0 * (l+1) * r0**l
        if verbose:
            print("pc = {:.8e} hc = {:.8e} dh = {:.8e} h0  = {:.8e}".format(pc,hc,dh,h0))
        return y, h0, h1
    
    def solve(self,pc,dh_fact=-1e-3):
        """
        Solves the Tolman-Oppenheimer-Volkov stellar structure equations using the pseudo-enthalpy formalism introduced in:
        Lindblom (1992) "Determining the Nuclear Equation of State from Neutron-Star Masses and Radii", Astrophys. J. 398 569.
        """
        # Initial data
        y, h0, h1 = self.__initial_data(pc, dh_fact=dh_fact)
        # Integrate
        sol = solve_ivp(self.__tov_rhs, [h0, h1], y,
                        max_step = h1, 
                        method = self.ode_method,
                        rtol = self.ode_rtol,
                        atol = self.ode_atol)
        # Take one final Euler step to get to surface 
        y = sol.y[:,-1]
        dy = self.__tov_rhs(sol.t[-1],y)
        y[:] -= dy[:] * h1
        sol.y.append(y)
        # Mass, Radius & Compactness
        M,R,C = self.__compute_mass_radius(y)
        # Match to Schwarzschild exterior
        sol.y[self.ivar['nu'],:] += np.log(1.0-(2.*M)/R) - sol.y[self.ivar['nu'],-1]
        # Even Love number k2 (for testing purposes only)
        #yyl = R*y[self.ivar['dH2']]/y[self.ivar['H2']]
        #k2 = self.__compute_Love_even_ell2(self, self.C,yyl)
        # Even Love numbers
        k = {}
        for l in self.leven:
            yyl = R*y[self.ivar['dH{}'.format(l)]]/y[self.ivar['H{}'.format(l)]]
            k[l] = self.__compute_Love_even(l,self.C,yyl)                   
        # Odd Love numbers
        j = {}
        for l in self.lodd:
            yyl = R*y[self.ivar['dPsi{}'.format(l)]]/y[self.ivar['Psi{}'.format(l)]]
            j[l] self.__compute_Love_odd(l,self.C,yyl)
        # Dump output
        if self.output:
            self.__output_solution(self,self.output,sol)
        return M,R,C,k,j

    def __compute_mass_radius(self, y):
        """
        Compute mass, radius, & compactness
        """
        R = y[self.ivar['r']]
        M = y[self.ivar['m']]
        return M,R,M/R

    def __compute_baryon_mass(self, sol):
        """
        Compute baryon mass
        """
        r = sol.y[self.ivar['r'],:]
        m = sol.y[self.ivar['m'],:]
        e = self.EOSEnergyDensityOfPseudoEnthalpyGeometerized(sol.t,eos)
        return np.trapz( 4*np.pi*r**2.*e./np.sqrt(1-2*m./r), r );

    def __compute_proper_radius(self, sol):
        """
        Compute baryon mass
        """
        r = sol.y[self.ivar['r'],:]
        m = sol.y[self.ivar['m'],:]
        return np.trapz( r, 1./np.sqrt((1-2*m./r)), r );
        
    def __compute_Love_even_ell2(self, c,y):
        """
        Eq. (50) of Damour & Nagar, Phys. Rev. D 80 084035 (2009).
        """
        a = 1. - 2. * c
        a2 = a**2
        c2 = c**2
        c3 = c * c2
        c4 = c * c3
        c5 = c * c4
        num = (8.0 / 5.0) * a2 * c5 * (2 * c * (y - 1) - y + 2)
        den = 2 * c * (4 * (y + 1) * c4 + (6 * y - 4) * c3 + (26 - 22 * y) * c2 + 3 * (5 * y - 8) * c - 3 * y + 6)
        den -= 3 * a2 * (2 * c * (y - 1) - y + 2) * np.log(1.0/a)
        return num / den;

    def __compute_Love_odd(self,ell,c,y):
        """
        Compute odd parity Love numbers given 
        * the multipolar index ell
        * the compactness c
        * the ratio y = R Psi(R)'/Psi(R) 
        Eq.(61) of Damour & Nagar, Phys. Rev. D 80 084035 (2009)
        """
        c2 = c**2
        c3 = c*c2
        c4 = c*c3
        c5 = c*c4
        j = 0.
        if ell == 2:
            nj =  96*c5*(-1 + 2*c)*(-3 + y)
            dj =  5.*(2*c*(9 + 3*c*(-3 + y) + 2*c2*(-3 + y) + 2*c3*(-3 + y) - 3*y + 12*c4*(1 + y)) + 3*(-1 + 2*c)*(-3 + y)*log(1 - 2*c))
            j = nj/dj
        return j
    
    def __compute_Love_even(self,ell,c,y):
        """
        Compute even parity Love numbers given 
        * the multipolar index ell
        * the compactness c
        * the ratio y = R H(R)'/H(R) 
        Eq.(49) of Damour & Nagar, Phys. Rev. D 80 084035 (2009)
        """
        c2 = c**2
        c3 = c*c2
        c4 = c*c3
        c5 = c*c4
        c6 = c*c5
        c7 = c*c6
        c8 = c*c7
        c9 = c*c8
        c10 = c*c9
        c11 = c*c10
        c13 = c2*c11
        c15 = c2*c13
        c17 = c2*c15
        k = 0.
        if ell < 2: return k
        if ell == 2:
            nk = (1-2*c)**2*(2+2*c*(y-1)-y)
            dk = 2*c*(6-3*y+3*c*(5*y-8))+4*c3*(13-11*y+c*(3*y-2)+2*c2*(1+y)) + 3*(1-2*c)**2*(2-y+2*c*(y-1))*np.log(1-2*c)
            k = 8*c5/5*nk/dk
        elif ell == 3:
            nk = (1 - 2*c)**2*(-3 - 3*c*(-2 + y) + 2*c2*(-1 + y) + y)
            dk = 2*c*(15*(-3 + y) + 4*c5*(1 + y) - 45*c*(-5 + 2*y) - 20*c3*(-9 + 7*y) + 2*c4*(-2 + 9*y) + 5*c2*(-72 + 37*y)) - 15*(1 - 2*c)**2*(-3 - 3*c*(-2 + y) + 2*c**2*(-1 + y) + y)*np.log(1.0/(1 - 2*c))
            k = 8*c7/7*nk/dk
        elif ell == 4:
            nk = (1 - 2*c)**2*(-7*(-4 + y) + 28*c*(-3 + y) - 34*c2*(-2 + y) + 12*c3*(-1 + y))
            dk = (2*c*(c2*(5360 - 1910*y) + c4*(1284 - 996*y) - 105*(-4 + y) + 8*c6*(1 + y) + 105*c*(-24 + 7*y) + 40*c3*(-116 + 55*y) + c5*(-8 + 68*y)) - 15*(1 - 2*c)**2*(-7*(-4 + y) + 28*c*(-3 + y) - 34*c2*(-2 + y) + 12*c3*(-1 + y))*np.log(1.0/(1 - 2*c)))
            k = 32*c9/147*nk/dk
        elif ell == 5:
            nk = (32*(1 - 2*c)**2*c11*(3*(-5 + y) - 15*c*(-4 + y) + 26*c2*(-3 + y) - 18*c3*(-2 + y) + 4*c4*(-1 + y)))
            dk = 99.*(2*c*(315*(-5 + y) + 8*c7*(1 + y) - 315*c*(-35 + 8*y) + 4*c6*(-2 + 27*y) - 56*c5*(-60 + 47*y) - 210*c3*(-170 + 57*y) + 105*c2*(-278 + 75*y) + 56*c4*(-345 + 158*y)) - 105*(1 - 2*c)**2*(3*(-5 + y) - 15*c*(-4 + y) + 26*c2*(-3 + y) - 18*c3*(-2 + y) + 4*c4*(-1 + y))*np.log(1.0/(1 - 2*c)))
            k = nk/dk
        elif ell == 6:
            nk = (1024*(1 - 2*c)**2*c13*(-33*(-6 + y) + 198*c*(-5 + y) - 444*c2*(-4 + y) + 456*c3*(-3 + y) - 208*c4*(-2 + y) + 32*c5*(-1 + y)))
            dk = 14157.*(2*c*(-3465*(-6 + y) + 32*c8*(1 + y) + 10395*c*(-16 + 3*y) + 16*c7*(-2 + 39*y) + 2016*c5*(-122 + 55*y) - 64*c6*(-457 + 362*y) - 210*c2*(-2505 + 541*y)  + 210*c3*(-3942 + 1015*y) - 84*c4*(-7917 + 2567*y)) - 105*(1 - 2*c)**2*(-33*(-6 + y) + 198*c*(-5 + y) - 444*c2*(-4 + y) + 456*c3*(-3 + y) - 208*c4*(-2 + y) + 32*c5*(-1 + y))*np.log(1.0/(1 - 2*c)))
            k = nk/dk
        elif ell == 7:
            nk = 1024*(1 - 2*c)**2*c15*(143*(-7 + y) - 1001*c*(-6 + y) + 2750*c2*(-5 + y) - 3740*c3*(-4 + y) + 2600*c4*(-3 + y) - 848*c5*(-2 + y) + 96*c6*(-1 + y))
            dk = 20449.*(2*c*(45045*(-7 + y) + 160*c9*(1 + y) - 45045*c*(-63 + 10*y) + 80*c8*(-2 + 53*y) - 432*c7*(-651 + 521*y) - 4620*c3*(-4333 + 902*y)  + 1155*c2*(-9028 + 1621*y) + 96*c6*(-33964 + 15203*y) + 126*c4*(-168858 + 42239*y) - 84*c5*(-144545 + 45971*y)) - 315*(1 - 2*c)**2*(143*(-7 + y) - 1001*c*(-6 + y) + 2750*c2*(-5 + y) - 3740*c3*(-4 + y) + 2600*c4*(-3 + y) - 848*c5*(-2 + y) + 96*c6*(-1 + y))*np.log(1.0/(1 - 2*c)))
            k = nk/dk
        elif ell == 8:
            nk = (16*(1 - 2*c)**2*(2*c*(c*(2*c*(2*c*(-737*(-4 + y) + 374*c*(-3 + y) - 92*c2*(-2 + y) + 8*c**3*(-1 + y)) + 1573*(-5 + y)) - 1859*(-6 + y)) + 572*(-7 + y)) - 143*(-8 + y)))
            dk = (286*c*(4*(90090 + c*(-900900 + c*(3768765 + c*(-8528520 + c*(11259633 + 2*c*(-4349499 + c*(1858341 + 8*c*(-47328 + c*(3092 + (-1 + c)*c))))))))) + (-1 + c)*(-1 + 2*c)*(-45045 + 4*c*(90090 + c*(-285285 + c*(450450 + c*(-365211 + 4*c*(34881 + c*(-4887 + 2*c*(36 + c))))))))*y) - 45045*(1 - 2*c)**2*(2*c*(c*(2*c*(2*c*(-737*(-4 + y) + 374*c*(-3 + y) - 92*c2*(-2 + y) + 8*c**3*(-1 + y)) + 1573*(-5 + y)) - 1859*(-6 + y)) + 572*(-7 + y)) - 143*(-8 + y))*np.log(1.0/(1 - 2*c)))
            k = 256/2431 * c17 * nk/dk
        else:
            #TODO k = ComputeLegendre(c,y,l) # https://bitbucket.org/bernuzzi/tov/src/master/ComputeLegendre.m
        return k
        
    def __compute_shape(self,ell,c,y,eps):
        """
        Compute even shape numbers given 
        * the multipolar index ell
        * the compactness c
        * the ratio y = R H(R)'/H(R) 
        Eq.(95) of Damour & Nagar, Phys. Rev. D 80 084035 (2009)
        """
        c2 = c**2
        c3 = c*c2
        c4 = c*c3
        c5 = c*c4
        c6 = c*c5
        c7 = c*c6
        c8 = c*c7
        c9 = c*c8
        c10 = c*c9
        c11 = c*c10
        c13 = c2*c11
        c15 = c2*c13
        c17 = c2*c15
        h = 0.
        if ell < 2: return h
        if ell == 2:
            nh = (-2 + 6*c + 2*c3*(1 + y) - c2*(6 + y))
            dh = (2*c*(6 + c2*(26 - 22*y) - 3*y + 4*c4*(1 + y) + 3*c*(-8 + 5*y) + c3*(-4 + 6*y)) - 3*(1 - 2*c)**2*(2 + 2*c*(-1 + y) - y)*np.log(1.0/(1 - 2*c)))
            h = -8*c5*nh/dh
        elif ell == 3:
            nh = -5 + 15*c + 2*c3*(1 + y) - c2*(12 + y)
            dh = (5.*(2*c*(15*(-3 + y) + 4*c5*(1 + y) - 45*c*(-5 + 2*y) - 20*c3*(-9 + 7*y) + 2*c4*(-2 + 9*y) + 5*c2*(-72 + 37*y)) - 15*(1 - 2*c)**2*(-3 - 3*c*(-2 + y) + 2*c2*(-1 + y) + y)*np.log(1.0/(1 - 2*c))))
            h = 16*c7*nh/dh
        elif ell == 4:
            nh = -9 + 27*c + 2*c3*(1 + y) - c2*(20 + y)
            dh = (21.*(2*c*(c2*(5360 - 1910*y) + c4*(1284 - 996*y) - 105*(-4 + y) + 8*c6*(1 + y) + 105*c*(-24 + 7*y)  + 40*c3*(-116 + 55*y) + c5*(-8 + 68*y)) - 15*(1 - 2*c)**2*(-7*(-4 + y) + 28*c*(-3 + y) - 34*c2*(-2 + y) + 12*c3*(-1 + y))*np.log(1.0/(1 - 2*c))))
            h = -64*c9*nh/dh
        return h
    
    def __compute_Love(self,ell,c,y,eps):
        """
        Compute Love numbers given 
        * the multipolar index ell
        * the compactness c
        * the ratio y = R H(R)'/H(R) or R Psi'(R)/Psi(R)
        * parity (0='even', 1 = 'odd')
        """
        if ell < 2: return 
        if eps % 2:
            # Odd
            return self.compute_Love_odd(ell,c,y)
        else:
            # Even
            return self.compute_Love_even(ell,c,y)

    def __compute_Lambda(self,ell,k,C):
        """
        Compute tidal polarizability $\Lambda_\ell$
        from Love numbers and compactness
        Note: Yagi's $\bar{\lambda}_\ell$ is $\Lambda_\ell$
        """
        div = 1.0/(factorial2(2*ell-1)*C**(2*ell+1))
        return 2.*kl*div
    
    def __output_makedir(self,outdir=None):
        """
        Make output dir
        """
        if not outdir: return
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        else:
            shutil.rmtree(outdir) # removes all the subdirs
            os.makedirs(outdir)
        return
        
    def __output_solution(self,outdir,sol):
        """
        Dump ODE solution
        """
        with open(outdir+"/tov.txt",'w') as f:
            np.savetxt(f,np.c_[sol.t,sol.y],fmt='%.9e',header=' '.join(self.var.key()), comments='#')
        return
    
    
    
    #TODO

    
    
    def __plot_mass_radius(self):
        """
        """
        return None

    def __output_sequence(self,outdir):
        """
        """
        return None
    
    def MakeSequence(self,pcrange):
        """
        Sequence from a range of values for the central pressure
        Optionally plot
        """
        return None
        
    def __find_target_global(self,key,val):
        """
        Conpute configuration with a target global parameter
        https://bitbucket.org/bernuzzi/tov/src/master/TOVIterate.m
        """
        return None
   

    # TOVSolver ---



if __name__ == "__main__":
    
    
    
    #TODO tests
    
    


    
