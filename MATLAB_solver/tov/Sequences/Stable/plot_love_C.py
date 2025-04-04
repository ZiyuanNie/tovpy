import matplotlib.pyplot as plt
import numpy as np
import glob

fig, ax = plt.subplots()
for f in glob.glob('*.txt'):

    C,k2 = np.loadtxt(f,usecols=(4,5),unpack=True)
    EOS = f.split('_')[0]
    ax.plot(C, k2,'o',label=EOS)

ax.set(xlabel=r'Compactness, $C=GM/(c^2R)$', ylabel=r'$k_2$',
       title='$\ell=2$ gravitoelectric Love number')
#ax.grid()
ax.legend(ncol=2)
fig.savefig("k2_C.png")
plt.show()
