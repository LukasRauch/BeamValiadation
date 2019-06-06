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


fig = plt.figure('Kragarm Benchmark',figsize=(6.6, 3.5))
plt.clf()
ax = fig.add_subplot(1, 1, 1)
major_ticks_x = np.arange(-1,1,0.1)
minor_ticks_x = np.arange(-1,1,0.05)
major_ticks_y = np.arange(-0,11,2)
minor_ticks_y = np.arange(-0,11,0.5)

ax.set_xticks(major_ticks_x)
ax.set_xticks(minor_ticks_x, minor=True)
ax.set_yticks(major_ticks_y)
ax.set_yticks(minor_ticks_y, minor=True)

# xx = np.array([1, 2, 3, 4, 5, 10, 15, 20])
# xx_lable = ['E1', 'E2', 'E3', 'E4', 'E5', 'E10', 'E15', 'E20']
xx = np.array([0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10])
w = np.array([[0.06636, 0.13098, 0.19235, 0.24945, 0.30172, 0.34901, 0.39147, 0.42941, 0.46326, 0.49346, 0.55566, 0.60325, 0.64039, 0.66996, 0.69397, 0.71379, 0.73042, 0.74457, 0.75676, 0.76737, 0.7767, 0.78498, 0.79239, 0.79906, 0.8051, 0.81061],
              [0.0663645, 0.1309754, 0.1923511, 0.2494538, 0.3017222, 0.3490122, 0.3914687, 0.4294168, 0.4632687, 0.4934632, 0.5556675, 0.6032638, 0.6404044, 0.669979, 0.6939841, 0.7138108, 0.7304458, 0.7445952, 0.756783, 0.7673972, 0.7767323, 0.7850149, 0.7924218, 0.7990925, 0.8051388, 0.8106505]])
u = np.array([[0.00265, 0.01035, 0.02249, 0.03817, 0.05643, 0.0764, 0.09732, 0.1186, 0.13981, 0.16064, 0.20996, 0.25442, 0.29394, 0.32894, 0.35999, 0.38763, 0.41236, 0.43459, 0.45468, 0.47293, 0.48957, 0.50483, 0.51886, 0.53182, 0.54383, 0.555],
             [0.0026465, 0.0103539, 0.0224878, 0.0381663, 0.0564336, 0.076401, 0.0973179, 0.1185978, 0.1398075, 0.1606443, 0.2099625, 0.2544256, 0.2939462, 0.3289491, 0.3599962, 0.3876383, 0.4123692, 0.4346013, 0.4546937, 0.4729406, 0.489588, 0.5048421, 0.518876, 0.531836, 0.5438461, 0.555012]])



labels = ["klassische Variation", "hyperduale Variation"]
# for i in range(2):
    # settings["color"] = colors[i]
plt.plot(w[0,:], xx, color='#005293', linestyle='-' , label = "$w^{Mattiasson}$")
plt.plot(w[1,:], xx, color='#005293', linestyle='none', marker='x',markersize=6,fillstyle='none', label = "$w^{IGA}$" )
plt.plot(u[0,:], xx, color='#A2AD00', linestyle='-' , label = "$u^{Mattiasson}$")
plt.plot(u[1,:], xx, color='#A2AD00', linestyle='none', marker='x',markersize=6,fillstyle='none', label = "$u^{IGA}$" )

# plt.plot(xx, dif, color='gray', label = 'Differenz')
plt.xlabel(r"Dimensionslose Verformung $u/L$, $w/L$")
plt.ylabel(r"Lastparameter $\symlambdaF")
# plt.xticks(xx, xx_lable)
plt.legend(loc='upper left')


# And a corresponding grid
# ax.grid(which='both')
ax.grid(which='minor', linestyle='dotted', alpha=0.4)
ax.grid(which='major', linestyle='dotted', alpha=0.8)

# plt.show()

tikz_save('D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/fig.tikz',
           figureheight = '\\figureheight',
           figurewidth = '\\figurewidth')

# # # plt.suptitle(r"$B$-spline basis functions, $p=%d$" % (3))
# plt.savefig("D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/kragarm_benchmark.pgf")
# plt.savefig("D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/kragarm_benchmark.pdf")

# plt.show()



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