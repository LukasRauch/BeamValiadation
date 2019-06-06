from __future__ import division, print_function, absolute_import
from pylab import *
from matplotlib2tikz import save as tikz_save

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors

# mpl.use("pgf")
# pgf_with_rc_fonts = {
#     "font.family": "serif",
#     "font.serif": [],                   # use latex default serif font
#     "font.sans-serif": ["DejaVu Sans"], # use a specific sans-serif font
#     "pgf.texsystem": "xelatex",
# }
# mpl.rcParams.update(pgf_with_rc_fonts)


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


# common settings for all plotted lines
settings = { "linestyle" : 'solid',
                "linewidth" : 1.5 }

# create a list of unique colors for plotting
#
# http://stackoverflow.com/questions/8389636/creating-over-20-unique-legend-colors-using-matplotlib
#
NUM_COLORS = nbasis = 3  # perform dummy evaluation to get number of basis functions
cm         = plt.get_cmap('gist_rainbow')
cNorm      = matplotlib.colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
scalarMap  = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
colors     = [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]


# fig = plt.figure('feder plot',figsize=(6.6, 3.5))
# plt.clf()
# ax = fig.add_subplot(1, 1, 1)
major_ticks_x = np.arange(0,1.1,0.2)
major_ticks_y = np.arange(0,1.3,0.2)
minor_ticks_x = np.arange(0,1.1,0.1)
minor_ticks_y = np.arange(-0,1.3,0.1)

# ax.set_xticks(major_ticks_x)
# ax.set_xticks(minor_ticks_x, minor=True)
# ax.set_yticks(major_ticks_y)
# ax.set_yticks(minor_ticks_y, minor=True)

markers = ['p','s','.','^','d','+','*',' ']
# markers = ['$a$','$b$','$c$','$d$','$e$','$f$','$g$','$h$']

# xx = np.array([1, 2, 3, 4, 5, 10, 15, 20])
# xx_lable = ['E1', 'E2', 'E3', 'E4', 'E5', 'E10', 'E15', 'E20']
xx = np.array([0, 0.027777778, 0.055555556, 0.083333333, 0.111111111, 0.138888889, 0.166666667, 0.194444444, 0.222222222, 0.25, 0.277777778, 0.305555556, 0.333333333, 0.361111111, 0.388888889, 0.416666667, 0.444444444, 0.472222222, 0.5, 0.527777778, 0.555555556, 0.583333333, 0.611111111, 0.638888889, 0.666666667, 0.694444444, 0.722222222, 0.75, 0.777777778, 0.805555556, 0.833333333, 0.861111111, 0.888888889, 0.916666667, 0.944444444, 0.972222222, 1])
xx_ref = np.array([0, 0.055555556, 0.111111111, 0.166666667, 0.222222222, 0.277777778, 0.333333333, 0.388888889, 0.444444444, 0.5, 0.555555556, 0.611111111, 0.666666667, 0.722222222, 0.777777778, 0.833333333, 0.888888889, 0.944444444, 1])

E1_x = np.array([0, 0.005069345, 0.020184943, 0.045071493, 0.079276793, 0.122181534, 0.173012675, 0.230860104, 0.294696046, 0.363397546, 0.435769596, 0.510570864, 0.586539581, 0.662419728, 0.736986779, 0.809072442, 0.877587819, 0.941544482, 1.000072995, 1.052438496, 1.098053011, 1.136484257, 1.167460791, 1.190873428, 1.206772963, 1.215364317, 1.216997285, 1.212154194, 1.201434799, 1.185538838, 1.16524671, 1.141398766, 1.114873757, 1.086566961, 1.057368534, 1.028142606, 0.999707611])
E1_y = np.array([0, 0.087045175, 0.172767994, 0.255872872, 0.335117163, 0.409336113, 0.477466056, 0.538565308, 0.59183221, 0.636620524, 0.672450279, 0.699015995, 0.716190459, 0.724024695, 0.722744021, 0.71274034, 0.694560877, 0.668893645, 0.636550036, 0.598444955, 0.555574998, 0.508995213, 0.459794983, 0.409073643, 0.357916361, 0.307370886, 0.258425668, 0.211989849, 0.168875571, 0.129782964, 0.095288137, 0.065834373, 0.041726704, 0.023129898, 0.010069842, 0.002438218, 2.74503E-07])
E1_z = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

E2_x = np.array([0, 0.005069345, 0.020184943, 0.045071493, 0.079276793, 0.122181534, 0.173012676, 0.230860105, 0.294696044, 0.363397547, 0.435769594, 0.510570869, 0.586539585, 0.662419728, 0.736986781, 0.809072447, 0.877587825, 0.941544488, 1.000072997, 1.052438495, 1.098053012, 1.136484259, 1.167460794, 1.190873428, 1.206772963, 1.215364316, 1.216997284, 1.212154193, 1.201434798, 1.185538837, 1.165246708, 1.141398763, 1.114873755, 1.086566959, 1.057368532, 1.028142604, 0.999707609])
E2_y = np.array([0, 0.087045175, 0.172767994, 0.255872873, 0.335117163, 0.409336113, 0.477466056, 0.538565308, 0.591832208, 0.636620524, 0.672450278, 0.699015996, 0.71619046, 0.724024695, 0.72274402, 0.712740339, 0.694560875, 0.668893643, 0.636550035, 0.598444955, 0.555574996, 0.50899521, 0.459794978, 0.409073641, 0.357916359, 0.307370884, 0.258425668, 0.211989847, 0.16887557, 0.129782963, 0.095288135, 0.065834372, 0.041726704, 0.023129898, 0.010069841, 0.002438218, 2.75E-07])
E2_z = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

E3_x = np.array([0, 0.005069345, 0.020184943, 0.045071493, 0.079276794, 0.122181535, 0.173012676, 0.230860105, 0.294696048, 0.363397558, 0.435769597, 0.510570863, 0.586539589, 0.662419726, 0.736986786, 0.809072463, 0.87758784, 0.941544491, 1.000072991, 1.052438498, 1.098053019, 1.136484259, 1.167460792, 1.190873429, 1.206772963, 1.215364316, 1.216997284, 1.212154193, 1.201434797, 1.185538836, 1.165246706, 1.141398763, 1.114873754, 1.086566959, 1.05736853, 1.028142603, 0.999707608])
E3_y = np.array([0, 0.087045175, 0.172767994, 0.255872873, 0.335117164, 0.409336114, 0.477466056, 0.538565308, 0.591832211, 0.636620531, 0.672450279, 0.699015994, 0.71619046, 0.724024694, 0.72274402, 0.712740336, 0.694560871, 0.668893641, 0.636550038, 0.598444952, 0.55557499, 0.50899521, 0.459794981, 0.409073637, 0.357916358, 0.307370885, 0.258425664, 0.211989845, 0.168875567, 0.129782961, 0.095288132, 0.06583437, 0.041726703, 0.023129897, 0.01006984, 0.002438218, 2.74487E-07])
E3_z = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])


ref_x = np.array([0, 0.020184464,0.079274571,0.173006657,0.294683402,0.435746721,0.586503328,0.736935592,0.877523058,1,1.097981554,1.167404623,1.206748336,1.217020492,1.201519028,1.165398669,1.115090679,1.057636208,1])
ref_y = np.array([0, 0.172767915, 0.33511665, 0.477464829, 0.591830898, 0.672450872, 0.716197244, 0.722763522, 0.694601255, 0.636619772, 0.555681004, 0.459940423, 0.358098622, 0.258634951, 0.169094542, 0.095492966, 0.041889581, 0.010162819, 0])
ref_z = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])


fig = plt.figure('Last-Verformungsauswertung Siralfeder')
ax = fig.add_subplot(1,1,1)
plt.tight_layout()

plt.xlabel('Lastfaktor $\lambda$')
plt.ylabel('Verschiebung: u, v, w')


plt.plot(xx_ref, ref_x, linewidth=0, linestyle = '',  marker="s", markersize=1,  color='black' , label='ref: u/L')
plt.plot(xx_ref, ref_y, linewidth=0, linestyle = '',  marker="s", markersize=1,  color='#005293' , label='ref: v/L')
plt.plot(xx_ref, ref_z, linewidth=0, linestyle = '',  marker="s", markersize=1,  color='#A2AD00' , label='ref: w/L')

plt.plot(xx, E1_x, linewidth=1.0, linestyle='--' , color='black' , label='#NP = 1: u/L')
plt.plot(xx, E1_y, linewidth=1.0, linestyle='--' , color='#005293' , label='#NP = 1: v/L')
plt.plot(xx, E1_z, linewidth=1.0, linestyle='--' , color='#A2AD00' , label='#NP = 1: w/L')


plt.plot(xx, E2_x, linewidth=1.0, linestyle='-.' , color='black' ,   label='#NP = 2: u/L')
plt.plot(xx, E2_y, linewidth=1.0, linestyle='-.' , color='#005293' , label='#NP = 2: v/L')
plt.plot(xx, E2_z, linewidth=1.0, linestyle='-.' , color='#A2AD00' , label='#NP = 2: w/L')


# plt.plot(xx, E2_x, linewidth=2.0, linestyle='-.' , color='#A2AD00' , label='E2 $\symDispu$')
# plt.plot(xx, E2_y, linewidth=2.0, linestyle='-.' , color='#A2AD00' , label='E2 $\symDispv$')
# plt.plot(xx, E2_z, linewidth=2.0, linestyle='-.' , color='#A2AD00' , label='E2 $\symDispw$')

plt.plot(xx, E3_x, linewidth=1.0, linestyle=':' , color='black' , label='#NP = 3: u/L')
plt.plot(xx, E3_y, linewidth=1.0, linestyle=':' , color='#005293' , label='#NP = 3: v/L')
plt.plot(xx, E3_z, linewidth=1.0, linestyle=':' , color='#A2AD00' , label='#NP = 3: w/L')

# plt.plot(xx, E3_x, linewidth=3.0, linestyle=':' , color='#005293' , label='E3 $\symDispu$')
# plt.plot(xx, E3_y, linewidth=3.0, linestyle=':' , color='#005293' , label='E3 $\symDispv$')
# plt.plot(xx, E3_z, linewidth=3.0, linestyle=':' , color='#005293' , label='E3 $\symDispw$')



plt.legend(bbox_to_anchor=(1, 1.0))
plt.grid(alpha=0.6, linestyle='dotted')
ax.set_xticks(major_ticks_x)
ax.set_yticks(major_ticks_y)



# # And a corresponding grid
# ax.grid(which='both')
# ax.grid(which='minor', linestyle='dotted', alpha=0.4)
# ax.grid(which='major', linestyle='dotted', alpha=0.8)


# tikz_save('D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/spiralfeder_plane.tikz',
#            figureheight = '\\figureheight',
#            figurewidth = '\\figurewidth')


# # # plt.suptitle(r"$B$-spline basis functions, $p=%d$" % (3))
# plt.savefig("D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/kragarm_benchmark.pgf")
# plt.savefig("D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/kragarm_benchmark.pdf")

plt.show()



# plt.figure('Polyn. Grad',figsize=(6.6, 3.5))
# plt.clf()

# xx_lable = ['D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10']
# xx = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10])
# yy = np.array([[0.032183421, 0.03218677, 0.032175384, 0.032657074, 0.033514649, 0.034861804, 0.036864568, 0.03966297, 0.043274208],
#               [0.002132842, 0.001909555, 0.001677701, 0.003348796, 0.004141974, 0.005118001, 0.005787625, 0.006215014, 0.007013678]])


# dif  = np.array([0.030050578, 0.030277215, 0.030497683, 0.029308278, 0.029372675, 0.029743802, 0.031076943, 0.033447956, 0.03626053])

# labels = ["klassische Variation", "hyperduale Variation"]
# for i in range(2):
#     settings["color"] = colors[i]
#     plt.plot( xx, yy[i,:], **settings , label = labels[i])
#     plt.plot( xx, yy[i,:], "kx" )

# plt.plot(xx, dif, color='gray', label='Differenz')
# plt.xlabel(u"Polynomialer Grad des NURBS-Patches")
# plt.ylabel(u"Berechnung der Steifigkeitsmatrix [s]")
# plt.xticks(xx, xx_lable)
# plt.legend(loc='center right')
# plt.grid(axis='y')
# # plt.suptitle(r"$B$-spline basis functions, $p=%d$" % (3))
# plt.savefig("D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/12_runtime_benchmark/poly_grad.pgf")
# plt.savefig("D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/12_runtime_benchmark/poly_grad.pdf")
# # plt.show()


# import pandas as pd
# data = pd.read_excel (r'D:\OneDrive\Uni\Masterarbeit\TimeMeasure\time.xlsx', sheet_name='poly_table') #for an earlier version of Excel, you may need to use the file extension of 'xls'
# df = pd.DataFrame(data, colums=[2,:])

# print (df)

print("done")