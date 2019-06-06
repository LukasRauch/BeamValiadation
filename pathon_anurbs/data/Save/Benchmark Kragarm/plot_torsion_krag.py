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
major_ticks_y = np.arange(0,180,60)
minor_ticks_y = np.arange(-0,11,0.5)

# ax.set_xticks(major_ticks_x)
# ax.set_xticks(minor_ticks_x, minor=True)
ax.set_yticks([0,60,120,180])
# ax.set_yticks(minor_ticks_y, minor=True)

# xx = np.array([1, 2, 3, 4, 5, 10, 15, 20])
# xx_lable = ['E1', 'E2', 'E3', 'E4', 'E5', 'E10', 'E15', 'E20']
xx = np.array([0, 0.125, 0.375, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.625, 3.875, 4, 4, 4.222222222, 4.666666667, 5.333333333, 5.777777778, 6, 6, 6.125, 6.375, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.625, 9.875, 10, ])

lam_02 = np.array([1.58E-06, 5.63E-01, 1.69E+00, 3.38E+00, 5.63E+00, 7.88E+00, 1.01E+01, 1.24E+01, 1.46E+01, 1.63E+01, 1.74E+01, 1.80E+01, 1.80E+01, 1.80E+01, 1.80E+01, 1.80E+01, 1.80E+01, 1.80E+01, 1.80E+01, 1.86E+01, 1.97E+01, 2.14E+01, 2.36E+01, 2.59E+01, 2.81E+01, 3.04E+01, 3.26E+01, 3.43E+01, 3.54E+01, 3.60E+01, ])
lam_04 = np.array([3.15E-06, 1.13E+00, 3.38E+00, 6.75E+00, 1.13E+01, 1.58E+01, 2.03E+01, 2.48E+01, 2.93E+01, 3.26E+01, 3.49E+01, 3.60E+01, 3.60E+01, 3.60E+01, 3.60E+01, 3.60E+01, 3.60E+01, 3.60E+01, 3.60E+01, 3.71E+01, 3.94E+01, 4.28E+01, 4.73E+01, 5.18E+01, 5.63E+01, 6.08E+01, 6.53E+01, 6.86E+01, 7.09E+01, 7.20E+01, ])
lam_06 = np.array([4.72E-06, 1.69E+00, 5.06E+00, 1.01E+01, 1.69E+01, 2.36E+01, 3.04E+01, 3.71E+01, 4.39E+01, 4.89E+01, 5.23E+01, 5.40E+01, 5.40E+01, 5.40E+01, 5.40E+01, 5.40E+01, 5.40E+01, 5.40E+01, 5.40E+01, 5.57E+01, 5.91E+01, 6.41E+01, 7.09E+01, 7.76E+01, 8.44E+01, 9.11E+01, 9.79E+01, 1.03E+02, 1.06E+02, 1.08E+02, ])
lam_08 = np.array([6.30E-06, 2.25E+00, 6.75E+00, 1.35E+01, 2.25E+01, 3.15E+01, 4.05E+01, 4.95E+01, 5.85E+01, 6.53E+01, 6.98E+01, 7.20E+01, 7.20E+01, 7.20E+01, 7.20E+01, 7.20E+01, 7.20E+01, 7.20E+01, 7.20E+01, 7.43E+01, 7.88E+01, 8.55E+01, 9.45E+01, 1.04E+02, 1.13E+02, 1.22E+02, 1.31E+02, 1.37E+02, 1.42E+02, 1.44E+02, ])
lam_10 = np.array([7.87E-06, 2.81E+00, 8.44E+00, 1.69E+01, 2.81E+01, 3.94E+01, 5.06E+01, 6.19E+01, 7.31E+01, 8.16E+01, 8.72E+01, 9.00E+01, 9.00E+01, 9.00E+01, 9.00E+01, 9.00E+01, 9.00E+01, 9.00E+01, 9.00E+01, 9.28E+01, 9.84E+01, 1.07E+02, 1.18E+02, 1.29E+02, 1.41E+02, 1.52E+02, 1.63E+02, 1.72E+02, 1.77E+02, 1.80E+02, ])


labels = ["klassische Variation", "hyperduale Variation"]
# for i in range(2):
    # settings["color"] = colors[i]
plt.plot(xx, lam_10,    color='black', linestyle='solid'  , label = "$\symlambdaF = 1.0$" )
plt.plot(xx, lam_08,    color='black', linestyle='dashed' , label = "$\symlambdaF = 0.8$" )
plt.plot(xx, lam_06,    color='black', linestyle='dashdot' , label = "$\symlambdaF = 0.6$")
plt.plot(xx, lam_04,    color='black', linestyle=(0, (3, 10, 1, 10, 1, 10))  , label = "$\symlambdaF = 0.4$" )
plt.plot(xx, lam_02,    color='black', linestyle='dotted'  , label = "$\symlambdaF = 0.2$")

# plt.plot(xx, dif, color='gray', label = 'Differenz')
plt.xlabel(r"Kragarm Laufkoordinate $\symx$ [m]")
plt.ylabel(r"Querschnittverdrehung $\symDispphi(\symx)$ [grad]")
# plt.xticks(xx, xx_lable)
plt.legend(loc='upper left')


# And a corresponding grid
# ax.grid(which='both')
ax.grid(linestyle='dotted', alpha=0.5)

# plt.show()

tikz_save('D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/09_Balkenimplementierung/tors_python.tikz',
           figureheight = '\\figureheight',
           figurewidth = '\\figurewidth')

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