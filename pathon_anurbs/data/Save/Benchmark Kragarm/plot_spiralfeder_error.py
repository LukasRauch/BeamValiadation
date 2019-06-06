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
major_ticks_y = np.arange(-0.0003,0.000051,0.00005)
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
xx = np.array([0, 0.027777778, 0.055555556, 0.083333333, 0.111111111, 0.138888889, 0.166666667, 0.194444444, 0.222222222, 0.25, 0.277777778, 0.305555556, 0.333333333, 0.361111111, 0.388888889, 0.416666667, 0.444444444, 0.472222222, 0.5, 0.527777778, 0.555555556, 0.583333333, 0.611111111, 0.638888889, 0.666666667, 0.694444444, 0.722222222, 0.75, 0.777777778, 0.805555556, 0.833333333, 0.861111111, 0.888888889, 0.916666667, 0.944444444, 0.972222222, 1])
# xx_ref = np.array([0, 0.055555556, 0.111111111, 0.166666667, 0.222222222, 0.277777778, 0.333333333, 0.388888889, 0.444444444, 0.5, 0.555555556, 0.611111111, 0.666666667, 0.722222222, 0.777777778, 0.833333333, 0.888888889, 0.944444444, 1])

E1_x = np.array([0, 1.14656E-07, 4.78869E-07, 1.15121E-06, 2.22221E-06, 3.80427E-06, 6.01837E-06, 8.97856E-06, 1.26444E-05, 1.73186E-05, 2.28752E-05, 2.92385E-05, 3.62529E-05, 4.36779E-05, 5.1187E-05, 5.83736E-05, 6.47613E-05, 6.98214E-05, 7.29947E-05, 7.37186E-05, 7.14573E-05, 6.57346E-05, 5.61677E-05, 4.24993E-05, 2.46275E-05, 2.63206E-06, -2.3207E-05, -5.23969E-05, -8.42293E-05, -0.000117787, -0.000151959, -0.000185472, -0.000216922, -0.000244826, -0.000267674, -0.000283988, -0.000292389])
E1_y = np.array([0, 1.16452E-08, 7.85631E-08, 2.41419E-07, 5.12787E-07, 8.65852E-07, 1.22623E-06, 1.46764E-06, 1.31204E-06, 7.51321E-07, -5.93039E-07, -3.0076E-06, -6.78483E-06, -1.22025E-05, -1.95011E-05, -2.88605E-05, -4.03774E-05, -5.40451E-05, -6.97359E-05, -8.71892E-05, -0.000106005, -0.000125645, -0.000145439, -0.000164604, -0.000182261, -0.000197475, -0.000209283, -0.000216742, -0.000218971, -0.000215204, -0.000204829, -0.000187441, -0.000162877, -0.00013125, -9.29769E-05, -4.87866E-05, 2.74503E-07])
E1_z = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

E2_x = np.array([0, 1.14668E-07, 4.78916E-07, 1.15131E-06, 2.22239E-06, 3.80455E-06, 6.01876E-06, 8.97907E-06, 1.26426E-05, 1.73192E-05, 2.28724E-05, 2.92438E-05, 3.62567E-05, 4.36784E-05, 5.11893E-05, 5.83786E-05, 6.47666E-05, 6.9827E-05, 7.29966E-05, 7.37181E-05, 7.14587E-05, 6.57359E-05, 5.61702E-05, 4.24994E-05, 2.46273E-05, 2.63163E-06, -2.32077E-05, -5.23979E-05, -8.42303E-05, -0.000117788, -0.000151961, -0.000185474, -0.000216924, -0.000244828, -0.000267676, -0.000283989, -0.000292391])
E2_y = np.array([0, 1.17304E-08, 7.87329E-08, 2.41673E-07, 5.13095E-07, 8.662E-07, 1.2266E-06, 1.46802E-06, 1.31045E-06, 7.51473E-07, -5.94498E-07, -3.00609E-06, -6.78437E-06, -1.22026E-05, -1.95014E-05, -2.88616E-05, -4.03793E-05, -5.40479E-05, -6.97374E-05, -8.71894E-05, -0.000106007, -0.000125647, -0.000145445, -0.000164606, -0.000182263, -0.000197477, -0.000209283, -0.000216743, -0.000218972, -0.000215205, -0.000204831, -0.000187443, -0.000162877, -0.000131251, -9.29773E-05, -4.87863E-05, 2.75E-07])
E2_z = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

E3_x = np.array([0, 1.14669E-07, 4.78951E-07, 1.15142E-06, 2.2226E-06, 3.80487E-06, 6.01915E-06, 8.97954E-06, 1.26469E-05, 1.73302E-05, 2.28755E-05, 2.92378E-05, 3.6261E-05, 4.36756E-05, 5.11948E-05, 5.83948E-05, 6.47818E-05, 6.98299E-05, 7.29911E-05, 7.37211E-05, 7.14653E-05, 6.57361E-05, 5.61683E-05, 4.25009E-05, 2.46275E-05, 2.63146E-06, -2.32077E-05, -5.23981E-05, -8.42309E-05, -0.000117789, -0.000151963, -0.000185475, -0.000216925, -0.000244828, -0.000267678, -0.000283991, -0.000292392])
E3_y = np.array([0, 1.17197E-08, 7.8869E-08, 2.41951E-07, 5.13504E-07, 8.66687E-07, 1.22705E-06, 1.46843E-06, 1.31375E-06, 7.58216E-07, -5.93172E-07, -3.00809E-06, -6.78361E-06, -1.2203E-05, -1.95017E-05, -2.88643E-05, -4.03839E-05, -5.40494E-05, -6.97346E-05, -8.71918E-05, -0.000106014, -0.000125648, -0.000145442, -0.000164609, -0.000182264, -0.000197476, -0.000209286, -0.000216745, -0.000218975, -0.000215207, -0.000204834, -0.000187444, -0.000162879, -0.000131251, -9.29783E-05, -4.8787E-05, 2.74487E-07])
E3_z = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])



fig = plt.figure('Last-Verformungsauswertung Siralfeder')
ax = fig.add_subplot(1,1,1)
plt.tight_layout()

plt.xlabel('Lastfaktor $\lambda$')
plt.ylabel('Fehler $u, v, w$')

plt.plot(xx, E1_x, linewidth=1.0, linestyle='--' , color='black'   , label='$\#$NP = 1: eps($\symDispu / \symLength$)')
plt.plot(xx, E1_y, linewidth=1.0, linestyle='--' , color='#005293' , label='$\#$NP = 1: eps($\symDispv / \symLength$)')
plt.plot(xx, E1_z, linewidth=1.0, linestyle='--' , color='#A2AD00' , label='$\#$NP = 1: eps($\symDispw / \symLength$)')


plt.plot(xx, E2_x, linewidth=1.0, linestyle='-.' , color='black'   , label='$\#$NP = 2: eps($\symDispu / \symLength$)')
plt.plot(xx, E2_y, linewidth=1.0, linestyle='-.' , color='#005293' , label='$\#$NP = 2: eps($\symDispv / \symLength$)')
plt.plot(xx, E2_z, linewidth=1.0, linestyle='-.' , color='#A2AD00' , label='$\#$NP = 2: eps($\symDispw / \symLength$)')

plt.plot(xx, E3_x, linewidth=1.5, linestyle=':' , color='black'   , label='$\#$NP = 3: eps($\symDispu / \symLength$)')
plt.plot(xx, E3_y, linewidth=1.5, linestyle=':' , color='#005293' , label='$\#$NP = 3: eps($\symDispv / \symLength$)')
plt.plot(xx, E3_z, linewidth=1.5, linestyle=':' , color='#A2AD00' , label='$\#$NP = 3: eps($\symDispw / \symLength$)')


plt.legend(loc='upper right',  prop={'size':10}, bbox_to_anchor=(1.8, 1))
plt.grid(alpha=0.5, linestyle='dotted')
ax.set_xticks(major_ticks_x)
ax.set_yticks(major_ticks_y)

# # And a corresponding grid
# ax.grid(which='both')
# ax.grid(which='minor', linestyle='dotted', alpha=0.4)
# ax.grid(which='major', linestyle='dotted', alpha=0.8)


# tikz_save('D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/spiralfeder_error.tikz',
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