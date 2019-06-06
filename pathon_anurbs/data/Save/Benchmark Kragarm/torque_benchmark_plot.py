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


# mpl.use("pgf")
# pgf_with_rc_fonts = {
#     "font.family": "serif",
#     "font.serif": [],                   # use latex default serif font
#     "font.sans-serif": ["DejaVu Sans"], # use a specific sans-serif font
#     'pgf.rcfonts': False,
# }
# mpl.rcParams.update(pgf_with_rc_fonts)

# plt.rcParams['svg.fonttype'] = 'none'



# import bspline
# import bspline.splinelab as splinelab


# # common settings for all plotted lines
# settings = { "linestyle" : 'solid',
#                 "linewidth" : 1.5 }

# # create a list of unique colors for plotting
# #
# # http://stackoverflow.com/questions/8389636/creating-over-20-unique-legend-colors-using-matplotlib
# #
NUM_COLORS = nbasis = 9  # perform dummy evaluation to get number of basis functions
cm         = plt.get_cmap('gist_rainbow')
cNorm      = matplotlib.colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
scalarMap  = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
colors     = [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]

def update_xlabels(ax):
    # xlabels = [format(label, '%d') for label in ax.get_xticks()]
    ax.set_xscale('log',basex=2)
    # ax.set_xticklabels(xlabels)

# fig = plt.figure('Kragarm Benchmark')
# plt.clf()

# f, ((ax1, ax2),(ax3,ax4)) = plt.subplots(2 , 2, sharex='row', sharey='none')
plt.figure('Torsionsbenchmark')
# update_xlabels(ax1)
# update_xlabels(ax2)
# update_xlabels(ax3)
# update_xlabels(ax4)
# major_ticks_x = np.arange(-1,1,0.1)
# minor_ticks_x = np.arange(-1,1,0.05)
# major_ticks_y = np.arange(-0,11,2)
# minor_ticks_y = np.arange(-0,11,0.5)

# ax.set_xticks(major_ticks_x)
# ax.set_xticks(minor_ticks_x, minor=True)
# ax.set_yticks(major_ticks_y)
# ax.set_yticks(minor_ticks_y, minor=True)

labels = ['$p = 2$', '$p = 3$', '$p = 4$', '$p = 5$', '$p = 6$', '$p = 7$', '$p = 8$', '$ref$']
xx_lable = ['1', '2', '4', '8', '16', '32', '64', '128', '256']
xx = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256])
disp = np.array([
[-0.037015628, -0.044570416, -0.081671702, -0.116652319, -0.132104011, -0.136730311, -0.137944419, -0.138251729, -0.1383288],
[-0.091011009, -0.129599369, -0.137672174, -0.138311406, -0.138351804, -0.138354331, -0.138354489, -0.138354479, -0.13835463],
[-0.137931492, -0.138283993, -0.138352375, -0.138354463, -0.138354499, -0.1383545, -0.138354505, -0.138354519, -0.138354517],
[-0.138271987, -0.13835251, -0.138354438, -0.138354498, -0.1383545, -0.138354501, -0.138354499, -0.138354452, -0.138354693],
[-0.138353953, -0.138354491, -0.1383545, -0.1383545, -0.1383545, -0.138354501, -0.138354508, -0.138354547, -0.138352705],
[-0.138354509, -0.1383545, -0.1383545, -0.1383545, -0.1383545, -0.138354499, -0.138354498, -0.138354429, -0.138356591],
[-0.138354503, -0.1383545, -0.1383545, -0.1383545, -0.1383545, -0.138354503, -0.138354474, -0.138354461, -0.138352782],
# reference
[-0.1383545, -0.1383545, -0.1383545, -0.1383545, -0.1383545, -0.1383545, -0.1383545, -0.1383545, -0.1383545]
]).transpose()

rotation = np.array([
[-0.002595279, 0.015413241, 0.041082709, 0.062874255, 0.072117581, 0.074846744, 0.075560365, 0.075740847, 0.0757861],
[0.11204205, 0.082439894, 0.076305706, 0.07582406, 0.07580151, 0.075801087, 0.075801174, 0.075801176, 0.075801293],
[0.076612478, 0.075869926, 0.075804552, 0.075801268, 0.075801204, 0.075801203, 0.075801201, 0.075801198, 0.075801032],
[0.075950072, 0.075814411, 0.075801136, 0.07580119, 0.07580121, 0.075801209, 0.075801189, 0.075801155, 0.075800843],
[0.075807183, 0.075801521, 0.075801182, 0.075801194, 0.075801201, 0.075801198, 0.075801212, 0.075801211, 0.075799907],
[0.07580076, 0.075801169, 0.075801196, 0.075801199, 0.075801206, 0.075801192, 0.075801209, 0.075801144, 0.075802891],
[0.075801264, 0.075801178, 0.075801206, 0.075801203, 0.075801209, 0.075801192, 0.075801184, 0.075801368, 0.075800023],
# reference
[0.0758012, 0.0758012, 0.0758012, 0.0758012, 0.0758012, 0.0758012, 0.0758012, 0.0758012, 0.0758012]
]).transpose()

bend = np.array([
[-0.003656655,  -0.005528608, -0.006890428,    -0.005741133,   -0.003816517,  -0.002488332, -0.001757968,    -0.00138155,    -0.001191311],
[-0.004481119,  -0.003320901, -0.001794533,    -0.001209523,   -0.001052584,  -0.001013078, -0.001003257,    -0.001000811,   -0.001000209],
[-0.00133283,   -0.000894676,  -0.000968105, -0.000995612,    -0.000999448,   -0.00099993,   -0.000999997,  -0.001000001, -0.001000011],
[-0.000741552,  -0.000941817, -0.000998597,    -0.000999947,   -0.000999994,  -0.001000001, -0.000999999,    -0.001000009,   -0.001000004],
[-0.000984201,  -0.00100051,  -0.001000067, -0.001000002,    -0.001, -0.001000004,    -0.000999993,   -0.001000013,  -0.000999989],
[-0.001002838,  -0.001000429, -0.000999961,    -0.000999995,   -0.000999993,  -0.000999999, -0.000999988,    -0.000999992,   -0.001000013],
[-0.001000165,  -0.00100001,  -0.000999994, -0.000999997,    -0.001000001,   -0.000999999,  -0.001000007, -0.001000004,    -0.000999993],
# referece
[-0.001, -0.001, -0.001, -0.001, -0.001, -0.001, -0.001, -0.001, -0.001]
]).transpose()


tors = np.array([
[-0.00019508,   -0.000306747,  -0.000583302, -0.000838073,    -0.000952592,   -0.000987555,  -0.000996841, -0.000999206,    -0.000999801],
[-0.000449666,  -0.000762327, -0.00096797, -0.000995674,    -0.000999447,   -0.000999931,  -0.000999991, -0.000999999,    -0.001000001],
[-0.000976698,  -0.001009885, -0.001001093,    -0.001000056,   -0.001000003,  -0.001,   -0.001,    -0.001, -0.001],
[-0.001021097,  -0.001003848, -0.000999852,    -0.000999995,   -0.001,    -0.001, -0.001,  -0.001,   -0.001000002],
[-0.001001012,  -0.000999965, -0.000999999,    -0.001, -0.001,  -0.001,   -0.001,    -0.001, -0.000999987],
[-0.000999834,  -0.000999979, -0.001000002,    -0.001, -0.001,  -0.001,   -0.001,    -0.001, -0.001000014],
[-0.000999992,  -0.001,   -0.001,    -0.001, -0.001,  -0.001,   -0.001,    -0.000999999,   -0.000999988],
# reference
[-0.001, -0.001, -0.001, -0.001, -0.001, -0.001, -0.001, -0.001, -0.001]
]).transpose()


markers = ['x','s','.','^','+','d','*',' ']

plt.subplot(2,2,1)
plt.xlabel('Anzahl Elemente', fontsize=14)
plt.ylabel('z-Verschiebung [m]', fontsize=14)
lineObjects = plt.plot(xx, disp, color='#005293', linewidth=1.0, markersize=5)
for i, curve in enumerate(lineObjects):
    curve.set_marker(markers[i])

lineObjects[-1].set_linewidth(2)
lineObjects[-1].set_color('#A2AD00')
plt.xscale('log',basex=2)
plt.xticks(xx, xx_lable)
plt.legend(iter(lineObjects), labels, loc='upper right ' ,  prop={'size':10})
# plt.ylim(-0.14, -0.12)
plt.grid(alpha=0.5)


plt.subplot(2,2,2)
plt.xlabel('Anzahl Elemente', fontsize=14)
plt.ylabel('$\phi$-Rotation [rad]', fontsize=14)
lineObjects = plt.plot(xx,rotation, color='#005293', linewidth=1.0, markersize=5, marker='s')
for i, curve in enumerate(lineObjects):
    curve.set_marker(markers[i])

lineObjects[-1].set_linewidth(2)
lineObjects[-1].set_color('#A2AD00')
plt.xscale('log',basex=2)
plt.xticks(xx, xx_lable)
plt.legend(iter(lineObjects), labels, loc='lower right' , prop={'size':10})
plt.grid(alpha=0.5)

plt.subplot(2,2,3)
plt.xlabel('Anzahl Elemente', fontsize=14)
plt.ylabel('Biegemoment [kNm]', fontsize=14)
lineObjects = plt.plot(xx,bend*1000, color='#005293', linewidth=1.0, markersize=5, marker='s')
for i, curve in enumerate(lineObjects):
    curve.set_marker(markers[i])

lineObjects[-1].set_linewidth(2)
lineObjects[-1].set_color('#A2AD00')
plt.xscale('log',basex=2)
plt.xticks(xx, xx_lable)
plt.legend(iter(lineObjects), labels, loc='lower right' , prop={'size':10})
plt.grid(alpha=0.5)

plt.subplot(2,2,4)
plt.xlabel('Anzahl Elemente', fontsize=14)
plt.ylabel('Torsionsmoment [kNm]', fontsize=14)
lineObjects = plt.plot(xx,tors*1000, color='#005293', linewidth=1.0, markersize=5, marker='s')
for i, curve in enumerate(lineObjects):
    curve.set_marker(markers[i])

lineObjects[-1].set_linewidth(2)
lineObjects[-1].set_color('#A2AD00')
plt.xscale('log',basex=2)
plt.xticks(xx, xx_lable)
plt.legend(iter(lineObjects), labels, loc='upper right' , prop={'size':10})
plt.grid(alpha=0.5)

plt.show()

# labels = ["klassische Variation", "hyperduale Variation"]
# # for i in range(2):
#     # settings["color"] = colors[i]
# plt.plot(w[0,:], xx, color='#005293', linestyle='-' , label = "$w^{Mattiasson}$")
# plt.plot(w[1,:], xx, color='#005293', marker='s',markersize=4,fillstyle='none', label = "$w^{IGA}$" )
# plt.plot(u[0,:], xx, color='#A2AD00', linestyle='-' , label = "$u^{Mattiasson}$")
# plt.plot(u[1,:], xx, color='#A2AD00', marker='s',markersize=4,fillstyle='none', label = "$u^{IGA}$" )

# # plt.plot(xx, dif, color='gray', label = 'Differenz')
# plt.xlabel(r"Dimensionslose Verformung $u/L$, $w/L$")
# plt.ylabel(r"Lastparameter $\symlambdaF")
# plt.legend(loc='upper left')


# # And a corresponding grid
# # ax.grid(which='both')
# ax.grid(which='minor', linestyle='dotted', alpha=0.4)
# ax.grid(which='major', linestyle='dotted', alpha=0.8)

plt.show()

# tikz_save('D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/torque_bench.tikz',
#            figureheight = '\\figureheight',
#            figurewidth = '\\figurewidth')

# # # plt.suptitle(r"$B$-spline basis functions, $p=%d$" % (3))
# plt.savefig("D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/torque_benchmark.pgf")
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