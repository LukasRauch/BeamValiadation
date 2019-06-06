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
# major_ticks_x = np.arange(-1,1,0.1)
# minor_ticks_x = np.arange(-1,1,0.05)
# major_ticks_y = np.arange(-0,11,2)
# minor_ticks_y = np.arange(-0,11,0.5)

# ax.set_xticks(major_ticks_x)
# ax.set_xticks(minor_ticks_x, minor=True)
# ax.set_yticks(major_ticks_y)
# ax.set_yticks(minor_ticks_y, minor=True)

markers = ['p','s','.','^','d','+','*',' ']
# markers = ['$a$','$b$','$c$','$d$','$e$','$f$','$g$','$h$']

# xx = np.array([1, 2, 3, 4, 5, 10, 15, 20])
# xx_lable = ['E1', 'E2', 'E3', 'E4', 'E5', 'E10', 'E15', 'E20']
xx = np.array([0.2, 0.4,  0.6,  0.8, 1.0])
c1_1 = np.array([0.199994, 0.399994, 0.599993, 0.799993, 0.999993])
c1_10 = np.array([0.0199665, 0.0399615, 0.0599604, 0.0799602, 0.0999601])
c1_1000 = np.array([0.00199873, 0.00399861, 0.00599859, 0.00799859, 0.00999859])
c1_10000 = np.array([0.000199888, 0.000399878, 0.000599876, 0.000799876, 0.000999876])


fig = plt.figure('Knickwinkel Kopplung')
ax = fig.add_subplot(1,1,1)

plt.xlabel('Lastfaktor $\lambda$F [MNm]')
plt.ylabel('Relativer Knickwinkel $\symPsihat$ [rad]')


plt.plot(xx, c1_1, linewidth=1.0,  marker=markers[0], markersize=5, color='black' , label='$C_k$ = 1e-0')
plt.plot(xx, c1_10, linewidth=1.0,  marker=markers[1], markersize=4, color='black' , label='$C_k$ = 1e-1')
plt.plot(xx, c1_1000, linewidth=1.0,  marker=markers[2], markersize=5, color='black' , label='$C_k$ = 1e-2')
plt.plot(xx, c1_10000, linewidth=1.0,  marker=markers[3], markersize=5, color='black' , label='$C_k$ = 1e-3')

plt.legend(loc='upper right',  prop={'size':10}, bbox_to_anchor=(1.5, 1))
plt.grid(alpha=0.5, linestyle='dotted')
ax.set_xscale('log')
ax.set_yscale('log')

# # And a corresponding grid
# # ax.grid(which='both')
# ax.grid(which='minor', linestyle='dotted', alpha=0.4)
# ax.grid(which='major', linestyle='dotted', alpha=0.8)


tikz_save('D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/09_Balkenimplementierung/krag_koppel.tikz',
           figureheight = '\\figureheight',
           figurewidth = '\\figurewidth/2-0.9')


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