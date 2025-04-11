from eos import *
from tovpy.units import Units
uts = Units()
from scipy.constants import c, G
import numpy as np
import os

EOS_FILE_NAME = [
    '2B.rns',
    '2H.rns',
    'ALF2.rns',
    'BHB_lp_25-Sept-2017.rns',
    'BLQ_30-Jan-2020.rns',
    'BLh_gibbs_180_0.35_10-Jul-2019.rns',
    'DD2_22-Jun-2018.rns',
    'DD2_25-Sept-2017.rns',
    'ENG.rns',
    'H4.rns',
    'LS220B0.rns',
    'LS220B0v2.rns',
    'LS_220_25-Sept-2017.rns',
    'MPA1.rns',
    'MS1.rns',
    'MS1b.rns',
    'NL3_05-Oct-2017.rns',
    'SFHo+BL_01-Apr-2019.rns',
    'SFHo_09-Feb-2019.rns',
    'SFHo_25-Sept-2017.rns',
    'SLy.rns',
    'SLy4.rns',
    'SLy4_15-Jun-2018.rns',
    'TM1_05-Oct-2017.rns',
    'TMA_05-Oct-2017.rns',
    'eosA',
    'eosAPR_fit',
    'eosAU',
    'eosB',
    'eosC',
    'eosF',
    'eosFPS_fit',
    'eosFPS_old',
    'eosG',
    'eosL',
    'eosN',
    'eosNV',
    'eosO',
    'eosSLy_fit',
    'eosUU',
    'eosWS',
    'eos_DD2_adb.rns',
    'eos_SFHo_adb.rns',
    'eos_SFHx_adb.rns'
]
# Note that LS1800B0.rns, eosAPR, eosFP, eosFPS, eosSLy, eosWNV are not included in the list above
from tov import TOV
pc = np.logspace(-12, -8, 200)
for eos_file in EOS_FILE_NAME:
    print(eos_file)
    M = np.zeros(len(pc))
    R = np.zeros(len(pc))
    C = np.zeros(len(pc))
    k2 = np.zeros(len(pc))
    j3 = np.zeros(len(pc))
    eos = EOS('tabular',name="from_file",filename=eos_file)
    this_tov = TOV(eos = eos,  leven = [2], 
                                         lodd = [3], 
                                         #ode_method='RK45',
                                         ode_atol=1e-10, 
                                         ode_rtol=1e-10, 
                                         dhfact=-1e-12
                                     )
    for i, pci in enumerate(pc):
        M[i],R[i],C[i],k,j = this_tov.solve(pci)
        print(pci)
        R[i] *= 1./1e3
        M[i] *= 1./uts.constant['MRSUN_SI'][0]
        k2[i] = k[2]
        j3[i] = j[3]
    data = np.column_stack((pc, M, R, C, k2, j3))
    module_dir = os.path.dirname(__file__)
    # save_path = os.path.join(module_dir, 'data', eos_file)
    save_path = os.path.join(module_dir, eos_file)
    np.savetxt(save_path + '.txt', data, delimiter='\t', 
           header='pc\tM\tR\tC\tk2\tj3', comments='')
    break
 