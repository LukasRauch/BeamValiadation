from __future__ import division, print_function, absolute_import
from pylab import *
from matplotlib2tikz import save as tikz_save

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors

mpl.use("pgf")
pgf_with_rc_fonts = {
    "font.family": "serif",
    "font.serif": [],                   # use latex default serif font
    "font.sans-serif": ["DejaVu Sans"], # use a specific sans-serif font
    'pgf.rcfonts': False,
}
mpl.rcParams.update(pgf_with_rc_fonts)

plt.rcParams['svg.fonttype'] = 'none'



import bspline
import bspline.splinelab as splinelab

major_ticks_x = np.arange(0,1.1,0.2)
major_ticks_y = np.arange(0,1.01,0.2)
# minor_ticks_x = np.arange(-1,1,0.05)
# minor_ticks_y = np.arange(-0,11,0.5)

# ax.set_xticks(major_ticks_x)
# ax.set_xticks(minor_ticks_x, minor=True)
# ax.set_yticks(major_ticks_y)
# ax.set_yticks(minor_ticks_y, minor=True)

markers = ['p','s','.','^','d','+','*',' ']
# markers = ['$a$','$b$','$c$','$d$','$e$','$f$','$g$','$h$']

# xx = np.array([1, 2, 3, 4, 5, 10, 15, 20])
# xx_lable = ['E1', 'E2', 'E3', 'E4', 'E5', 'E10', 'E15', 'E20']
xx = np.array([0, 0.055555556, 0.111111111, 0.166666667, 0.222222222, 0.277777778, 0.333333333, 0.388888889, 0.444444444, 0.5, 0.555555556, 0.611111111, 0.666666667, 0.722222222, 0.777777778, 0.833333333, 0.888888889, 0.944444444, 1])
xx_ref = np.array([0, 0.055555556, 0.111111111, 0.166666667, 0.222222222, 0.277777778, 0.333333333, 0.388888889, 0.444444444, 0.5, 0.555555556, 0.611111111, 0.666666667, 0.722222222, 0.777777778, 0.833333333, 0.888888889, 0.944444444, 1])

E1_x = np.array([0, 0.019086952, 0.066058656, 0.125507085, 0.189028841, 0.253166661, 0.316588379, 0.378814814, 0.439688723, 0.499301518, 0.557702143, 0.61501255, 0.671370185, 0.726924621, 0.781840973, 0.836309397, 0.890563168, 0.945000349, 1.000001702])
E1_y = np.array([0, -0.118843906, -0.21737093, -0.292984796, -0.350628388, -0.395026876, -0.429483566, -0.456203678, -0.476663315, -0.491903355, -0.502612169, -0.509272306, -0.512197584, -0.511561978, -0.507408623, -0.499638526, -0.487970891, -0.471822694, -0.450157014])
E1_z = np.array([0, 0.118843906, 0.21737093, 0.292984796, 0.350628388, 0.395026876, 0.429483566, 0.456203678, 0.476663315, 0.491903355, 0.502612169, 0.509272306, 0.512197584, 0.511561978, 0.507408623, 0.499638526, 0.487970891, 0.471822694, 0.450156922])

ref_x = np.array([0, 0.005069345, 0.020184943, 0.045071493, 0.079276793, 0.122181534, 0.173012675, 0.230860104, 0.294696046, 0.363397546, 0.435769596, 0.510570864, 0.586539581, 0.662419728, 0.736986779, 0.809072442, 0.877587819, 0.941544482, 1.000072995])
# ref_y = np.array([0, 0.172767915, 0.33511665, 0.477464829, 0.591830898, 0.672450872, 0.716197244, 0.722763522, 0.694601255, 0.636619772, 0.555681004, 0.459940423, 0.358098622, 0.258634951, 0.169094542, 0.095492966, 0.041889581, 0.010162819, 0])
# ref_z = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])


fig = plt.figure('Last-Verformungsauswertung Spiralfeder')
ax = fig.add_subplot(1,1,1)
plt.tight_layout()

plt.xlabel('Lastfaktor $\lambda$')
plt.ylabel('Verschiebung $u, v, w$')


# plt.plot(xx_ref, ref_x, linewidth=0, linestyle = '',  marker="s", markersize=1,  color='black' , label='ref: $\symDispu$/L')
# plt.plot(xx_ref, ref_y, linewidth=0, linestyle = '',  marker="s", markersize=1,  color='#005293' , label='ref: $\symDispv$/L')
# plt.plot(xx_ref, ref_z, linewidth=0, linestyle = '',  marker="s", markersize=1,  color='#A2AD00' , label='ref: $\symDispw$/L')

plt.plot(xx, E1_x, linewidth=1.0, linestyle='-' , color='black' , label='$\#$NP = 1: $\symDispu / \symLength$')
plt.plot(xx, E1_z, linewidth=2.5, linestyle='-' , color='#A2AD00' , label='$\#$NP = 1: $\symDispw / \symLength$')
plt.plot(xx, -E1_y, linewidth=2.0, linestyle=':' , color='#005293' , label='$\#$NP = 1: $-\symDispv / \symLength$')

# plt.plot(xx, ref_x, linewidth=1.0, linestyle='-.' , color='gray' , label='$\symDispu_{ref}$')


plt.legend(loc='upper right',  prop={'size':10}, bbox_to_anchor=(1.62, 1))
plt.grid(alpha=0.5, linestyle='dotted')
ax.set_xticks(major_ticks_x)
ax.set_yticks(major_ticks_y)

# # And a corresponding grid
# ax.grid(which='both')
# ax.grid(which='minor', linestyle='dotted', alpha=0.4)
# ax.grid(which='major', linestyle='dotted', alpha=0.8)


tikz_save('D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/spiralfeder_zwei_plane.tikz',
           figureheight = '\\figureheight',
           figurewidth = '\\figurewidth')


# # # plt.suptitle(r"$B$-spline basis functions, $p=%d$" % (3))
# plt.savefig("D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/kragarm_benchmark.pgf")
# plt.savefig("D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/kragarm_benchmark.pdf")

plt.show()

print("done")