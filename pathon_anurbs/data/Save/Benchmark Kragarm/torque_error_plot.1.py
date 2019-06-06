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

xx1 = np.array([1])
xx2 = np.array([1, 2])
xx3 = np.array([1, 2, 4])
xx4 = np.array([1, 2, 4, 8])
xx5 = np.array([1, 2, 4, 8, 16])
xx6 = np.array([1, 2, 4, 8, 16, 32])
xx7 = np.array([1, 2, 4, 8, 16, 32, 64])
xx8 = np.array([1, 2, 4, 8, 16, 32, 64, 128])
xx9 = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256])

xxr1 = np.array([256])
xxr2 = np.array([128, 256])
xxr3 = np.array([64, 128, 256])
xxr4 = np.array([32, 64, 128, 256])
xxr5 = np.array([16, 32, 64, 128, 256])
xxr6 = np.array([8, 16, 32, 64, 128, 256])
xxr7 = np.array([4, 8, 16, 32, 64, 128, 256])


disp = np.array([
[0.732458085,   0.677853514,   0.409692478,   0.156859233,   0.045177342,   0.01173933,    0.002963991,    0.000742812,    0.000185756],
[0.342189745,   0.063280421,   0.004931719,   0.000311473,   1.94833E-05,   1.21962E-06,   7.92731E-08],
[0.003057421,   0.000509608,   1.53567E-05,   2.68942E-07,   4.59593E-09,   1.06138E-09],
[0.00059639,    1.43869E-05,    4.50114E-07,    1.13292E-08,    5.397E-11],
[3.95121E-06,   6.65652E-08,   3.52934E-10,   2.79984E-10],
[6.12469E-08,   2.04119E-09,   2.63649E-10],
[2.16908E-08,   1.97348E-10,	6.87437E-11]
]).transpose()

disp_g = np.array([
[0],
[7.92731E-08,   1.51125E-07,   9.37278E-07],
[1.06138E-09,   3.67796E-08,   1.40336E-07,   1.23254E-07],
[5.397E-11,     7.96044E-09,   8.25041E-09,   3.44805E-07,   1.39153E-06],
[2.79984E-10,   9.59195E-10,   6.224E-09, 5.45298E-08, 3.35733E-07, 1.29715E-05],
[2.63649E-10,   2.6834E-10,    2.33839E-09,    9.2616E-09, 1.77353E-08, 5.14175E-07, 1.51111E-05],
[6.87437E-11,   5.49314E-11,   6.23514E-10,   1.89525E-08,   1.87923E-07,   2.7938E-07,    1.24194E-05]
]).transpose()

rotation = np.array([
[1.034237966,   0.796662308,   0.458020343,   0.170537481,   0.04859579,    0.012591567,    0.003177194,    0.000796205,    0.000199208],
[0.478103895,   0.087580327,   0.006655642,   0.000301578,   4.08053E-06,   1.49006E-06,   3.47956E-07,   3.23197E-07],
[0.010702695,   0.00090666,    4.4215E-05, 8.9839E-07,  5.35024E-08,  4.10551E-08,  5.64495E-09],
[0.00196398,    0.000174279,    8.44158E-07,    1.35654E-07,    1.2775E-07, 1.20269E-07],
[7.89287E-05,   4.23459E-06,   2.41989E-07,   8.32898E-08,   5.63687E-09],
[5.80502E-06,   4.13409E-07,   5.14538E-08,   2.25924E-08],
[8.43095E-07,   2.90728E-07,   7.05803E-08,   3.6302E-08]
]).transpose()

rotation_g = np.array([
[0],
[3.23197E-07,   1.22409E-06],
[5.64495E-09,   2.48929E-08,   2.21728E-06],
[1.20269E-07,   1.52012E-07,   5.95588E-07,   4.71335E-06],
[5.63687E-09,   3.14514E-08,   1.54974E-07,   1.34291E-07,   1.70613E-05],
[2.25924E-08,   7.33856E-08,   1.10337E-07,   1.09579E-07,   7.37348E-07,   2.23006E-05],
[3.6302E-08,    1.08193E-07,   1.13359E-07,   2.2173E-07,    2.21058E-06,    1.55338E-05]
]).transpose()

bend = np.array([
[2.656655153,   4.528607549,   5.890428459,   4.741133453,   2.816517378,   1.488331905,   0.757968483,   0.381550088,   0.191311184],
[3.48111939,    2.32090088, 0.794533071, 0.209523064, 0.052583729, 0.013077685, 0.003257248, 0.000810645, 0.000209211],
[0.332830312,   0.105324073,   0.031895053,   0.00438767,    0.000552181,    6.97504E-05,    2.74165E-06,    1.11927E-06],
[0.258448321,   0.058182576,   0.001402751,   5.3063E-05,    5.62481E-06,    1.29491E-06,    8.27488E-07],
[0.015798681,   0.000509802,   6.72766E-05,   1.95053E-06,   4.68008E-07],
[0.002837656,   0.000429069,   3.85199E-05,   5.36353E-06,   6.71235E-06,   6.33997E-07],
[0.000164892,   1.03934E-05,   5.55113E-06,   3.4763E-06,    1.37008E-06,    1.38575E-06]
]).transpose()

bend_g = np.array([
[0],
[0],
[1.11927E-06,   1.09958E-05],
[8.27488E-07,   8.75194E-06,   4.15072E-06],
[4.68008E-07,   3.58057E-06,   6.93975E-06,   1.3245E-05,    1.10106E-05],
[6.33997E-07,   1.23969E-05,   8.27455E-06,   1.32825E-05],
[1.38575E-06,   7.20416E-06,   4.33605E-06,   6.92791E-06]
]).transpose()


tors = np.array([
[0.804919983,0.69325323, 0.41669785,0.161927281, 0.047407992,0.012445367, 0.003159418, 0.000793977, 0.000198783],
[0.550333613,0.23767281, 0.032029672, 0.004326361, 0.000552828, 6.93572E-05, 8.67494E-06, 1.23302E-06, 6.74811E-07],
[0.023302188,   0.009884506,   0.001092974,   5.55813E-05,   2.63365E-06,   1.49861E-07,   2.8209E-08],
[0.021096801,   0.00384847,    0.000147958,    5.29668E-06,    1.5667E-07, 1.4E-10],
[0.001011881,   3.51579E-05,   7.18454E-07,   2.13071E-08],
[0.000165591,   2.10374E-05,   2.43568E-06,   6.3318E-08,    1.10813E-07,    2.4974E-08, 3.1222E-09],
[7.5633E-06,    4.25845E-07,    2.01538E-07,    1.04149E-07,    1.13553E-07,    3.8916E-08]
]).transpose()

tors_g = np.array([
[0],
[0],
[2.8209E-08,    1.74871E-07,   1.58754E-07],
[1.4E-10,       1.23816E-08,   3.11418E-07,   2.3653E-06],
[2.13071E-08,   1.05467E-07,   1.56337E-08,   5.8818E-08,    3.87744E-07,    1.26218E-05],
[3.1222E-09,    4.65967E-07,   1.36437E-05],
[3.8916E-08,    2.23879E-07,   7.8464E-07,    1.19847E-05]
]).transpose()


markers = ['x','s','.','^','+','d','*',' ']
lineObjects = []

plt.subplot(2,2,1)
# plt.hold()
plt.xlabel('Anzahl Elemente', fontsize=14)
plt.ylabel('$\mu_u$ z-Verschiebung', fontsize=14)

# gray
plt.plot(xxr3, disp_g[1], marker=markers[1], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr4, disp_g[2], marker=markers[2], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr5, disp_g[3], marker=markers[3], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr6, disp_g[4], marker=markers[4], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr7, disp_g[5], marker=markers[5], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr7, disp_g[6], marker=markers[6], linewidth=1.0, markersize=5, color='#CCCCC6')

plt.plot(xx9, disp[0], label='$p=2$', marker=markers[0], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx7, disp[1], label='$p=3$', marker=markers[1], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx6, disp[2], label='$p=4$', marker=markers[2], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx5, disp[3], label='$p=5$', marker=markers[3], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx4, disp[4], label='$p=6$', marker=markers[4], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx3, disp[5], label='$p=7$', marker=markers[5], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx3, disp[6], label='$p=8$', marker=markers[6], linewidth=1.0, markersize=5, color='#005293')

plt.xscale('log',basex=2)
plt.yscale('log',basey=10)
plt.xticks(xx9, xx_lable)
plt.legend(loc='upper right' ,  prop={'size':10})
# plt.ylim(-0.14, -0.12)
plt.grid(alpha=0.5)


plt.subplot(2,2,2)
plt.xlabel('Anzahl Elemente', fontsize=14)
plt.ylabel('$\mu_{\phi}$ Rotation', fontsize=14)

# gray
plt.plot(xxr2, rotation_g[1], marker=markers[1], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr3, rotation_g[2], marker=markers[2], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr4, rotation_g[3], marker=markers[3], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr5, rotation_g[4], marker=markers[4], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr6, rotation_g[5], marker=markers[5], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr6, rotation_g[6], marker=markers[6], linewidth=1.0, markersize=5, color='#CCCCC6')

plt.plot(xx9, rotation[0], label='$p=2$', marker=markers[0], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx8, rotation[1], label='$p=3$', marker=markers[1], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx7, rotation[2], label='$p=4$', marker=markers[2], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx6, rotation[3], label='$p=5$', marker=markers[3], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx5, rotation[4], label='$p=6$', marker=markers[4], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx4, rotation[5], label='$p=7$', marker=markers[5], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx4, rotation[6], label='$p=8$', marker=markers[6], linewidth=1.0, markersize=5, color='#005293')

plt.xscale('log',basex=2)
plt.yscale('log',basey=10)
plt.xticks(xx9, xx_lable)
plt.legend(loc='upper right' , prop={'size':10})
plt.grid(alpha=0.5)

plt.subplot(2,2,3)
plt.xlabel('Anzahl Elemente', fontsize=14)
plt.ylabel('$\mu_{M_y}$ Biegemoment', fontsize=14)

# gray
plt.plot(xxr2, bend_g[2], marker=markers[2], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr3, bend_g[3], marker=markers[3], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr5, bend_g[4], marker=markers[4], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr4, bend_g[5], marker=markers[5], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr4, bend_g[6], marker=markers[6], linewidth=1.0, markersize=5, color='#CCCCC6')

plt.plot(xx9, bend[0], label='$p=2$', marker=markers[0], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx9, bend[1], label='$p=3$', marker=markers[1], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx8, bend[2], label='$p=4$', marker=markers[2], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx7, bend[3], label='$p=5$', marker=markers[3], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx5, bend[4], label='$p=6$', marker=markers[4], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx6, bend[5], label='$p=7$', marker=markers[5], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx6, bend[6], label='$p=8$', marker=markers[6], linewidth=1.0, markersize=5, color='#005293')

plt.xscale('log',basex=2)
plt.yscale('log',basey=10)
plt.xticks(xx9, xx_lable)
plt.legend(loc='upper right' , prop={'size':10})
plt.grid(alpha=0.5)

plt.subplot(2,2,4)
plt.xlabel('Anzahl Elemente', fontsize=14)
plt.ylabel('$\mu_{M_T}$ Torsionsmoment', fontsize=14)

# gray
plt.plot(xxr3, tors_g[2], marker=markers[2], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr4, tors_g[3], marker=markers[3], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr6, tors_g[4], marker=markers[4], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr3, tors_g[5], marker=markers[5], linewidth=1.0, markersize=5, color='#CCCCC6')
plt.plot(xxr4, tors_g[6], marker=markers[6], linewidth=1.0, markersize=5, color='#CCCCC6')

plt.plot(xx9, tors[0], label='$p=2$', marker=markers[0], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx9, tors[1], label='$p=3$', marker=markers[1], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx7, tors[2], label='$p=4$', marker=markers[2], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx6, tors[3], label='$p=5$', marker=markers[3], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx4, tors[4], label='$p=6$', marker=markers[4], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx7, tors[5], label='$p=7$', marker=markers[5], linewidth=1.0, markersize=5, color='#005293')
plt.plot(xx6, tors[6], label='$p=8$', marker=markers[6], linewidth=1.0, markersize=5, color='#005293')
plt.xscale('log',basex=2)
plt.yscale('log',basey=10)
plt.xticks(xx9, xx_lable)
plt.legend(loc='upper right' , prop={'size':10})
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
# plt.savefig("D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/torque_result_TUM.pdf")

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