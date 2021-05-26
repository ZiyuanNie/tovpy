#!/usr/bin/env python3

class Units(object):

    """

    Class for physical constants and conversion factors
    https://lscsoft.docs.ligo.org/lalsuite/lal/_l_a_l_constants_8h_source.html

    """

    def __init__(self):

        # Physical constants
        
        self.constant = {
            'C_SI': (299792458e0, 'Speed of light in vacuum, m s^-1'),
            'H_SI': (6.62607015e-34, 'Planck constant, J s'),
            'QE_SI': (1.602176634e-19, 'Electron charge, C'), 
            'MOL': (6.02214076e+23, 'Avogadro constant, dimless'),
            'K_SI': (1.380649e-23, 'Boltzmann constant, J K^-1'),
            'R_SI': (8.31446261815324, 'Molar gas constant, J mol^-1 K^-1'), 
            'GEARTH_SI': (9.80665, 'Standard gravity, m s^-2'), 
            'PATM_SI': (101325e0, 'Standard atmosphere, Pa'), 
            'HBAR_SI': (1.054571817646156391262428003302280745e-34, 'Reduced Planck constant, J s'),
            'ALPHA': (0.0072973525693, 'Fine structure constant, dimless'),
            'RYD_SI': (10973731.568160, 'Rydberg constant, m^-1'),
            'MP_ME': (1836.15267343, 'Proton-electron mass ratio, dimless'),
            'ME_AMU': (0.000548579909065, 'Electron mass, atomic mass units'),
            'G_SI': (6.67430e-11, 'Gravitational constant, N m^2 kg^-2'),
            'AMU_SI': (1.660539066595378801332508797951914123e-27, 'Atomic mass unit, kg'),
            'MP_SI': (1.672621923684144692109494784075478798e-27, 'Proton mass, kg'),
            'ME_SI': (9.109383701517728819842163772087735080e-31, 'Electron mass, kg'),
            'MSUN_SI': (1.988409902147041637325262574352366540e30, 'Solar mass, kg'),
            'GMSUN_SI': (1.32712442099e+20, 'Solar mass parameter, m^3 s^-2 (TCB)'),
            'MRSUN_SI': (1.476625061404649406193430731479084713e3, 'Geometrized solar mass, m'),
            'MTSUN_SI': (4.925491025543575903411922162094833998e-6, 'Geometrized solar mass, s'),
            'NUCL_DENS_SI': (2.8e17, 'Nuclear density in kg m^-3'),
        }
        
        G_SI = self.constant['G_SI'][0]
        C_SI = self.constant['C_SI'][0]
        MSUN_SI = self.constant['MSUN_SI'][0]
        MRSUN_SI = self.constant['MRSUN_SI'][0]
        MTSUN_SI = self.constant['MTSUN_SI'][0]
        
        G_C2_SI = G_SI / C_SI**2 
        G_C4_SI = G_C2_SI / C_SI**2 
        NUCL_DENSI_GEOM_SI  = self.constant['NUCL_DENS_SI'][0] * G_C2_SI 

        self.constant['C_CGS'] = (self.constant['C_SI'][0] * 100, 'Speed of light in vacuum, cm s^-1')
        self.constant['G_C2_SI'] = (G_C2_SI, 'Geometrized density, kg/m^3')
        self.constant['G_C4_SI'] = (G_C4_SI, 'Geometrized pressure in Pa')
        self.constant['NUCL_DENSI_GEOM_SI'] = (NUCL_DENSI_GEOM_SI, 'Nuclear density in geometrized units of m^-2')
        
        # Some conversion factors from geometrized to CGS or SI unit systems
        #TODO check/fixme
        
        self.conversion_factor = {
            'pressure': {'cgs': 1.0 / G_C4_SI * 10., 'si': 1.0 / G_C4_SI, 'geom': 1.0},
            'energy_density': {'cgs': 1.0 / G_C4_SI * 10., 'si': C_SI ** 4. / G_SI, 'geom': 1.},
            'density': {'cgs': 1.0 / G_C2_SI / 1000., 'si': 1.0 / G_C2_SI , 'geom': 1.0},
            'pseudo_enthalpy': {'cgs': 1.0, 'si': 1.0, 'geom': 1.0},
            'mass': {'cgs': 1.0 / G_C2_SI * 1000., 'si': 1.0 / G_C2_SI, 'geom': 1.0, 'm_sol': 1.0 / G_C2_SI / MSUN_SI},
            'radius': {'cgs': MRSUN_SI * 100., 'si': MRSUN_SI, 'geom': 1.0, 'm_sol': 1.0 / MSUN_SI },
            'time': {'cgs': MTSUN_SI * 100., 'si': MTSUN_SI, 'geom': 1.0},
        }
        
    def show(self):
        """
        Screen the physical constants
        """
        for k,v in self.constant.items():
            print("{:<20} {:<50} {}".format(k,v[1],v[0]))
        return

    def const(self,k):
        """
        Return value of a constant given its keyword
        """
        return self.constant[k][0]

    def geom_to_gcs(self,k):
        """
        Return conversion factor from geom to CGS of a quantity k
        """
        return self.conversion_factor[k]['cgs']

    def geom_to_si(self,k):
        """
        Return conversion factor from geom to SI of a quantity k
        """
        return self.conversion_factor[k]['si']


if __name__ == "__main__":

    
    uts = Units()
    uts.show()

    #print(uts.const('MSUN_SI'))
    #print(uts.geom_to_si('mass'))
