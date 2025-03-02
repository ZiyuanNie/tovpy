from eos import EOS
EoS = EOS('tabular',name="from_file",filename='eosG')
from tov import TOV
tov = TOV(eos=EoS)
M,R,C,k,j = tov.solve(1e-9)[:-1]
print(M,R,C,k,j)
print(1-2*M/R)
# print(EoS.EnergyDensity_Of_Pressure(1e-9))
# print(not None)