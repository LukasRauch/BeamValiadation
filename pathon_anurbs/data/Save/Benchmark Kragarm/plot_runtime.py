from __future__ import division, print_function, absolute_import
from pylab import *
from matplotlib2tikz import save as tikz_save

import numpy as np
import matplotlib2tikz
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors

#Direct input 
# plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
# params = {'text.usetex' : True,
#           'font.size' : 11,
#           'font.family' : 'lmodern',
#           'text.latex.unicode': True,
#           }
# plt.rcParams.update(params) 



import bspline
import bspline.splinelab as splinelab

# fig = plt.figure('feder plot',figsize=(6.6, 3.5))
# plt.clf()
# ax = fig.add_subplot(1, 1, 1)
# major_ticks_x = np.arange(0,1.1,0.2)
# major_ticks_y = np.arange(0,1.3,0.2)
# # minor_ticks_x = np.arange(-1,1,0.05)
# minor_ticks_y = np.arange(-0,11,0.5)

# ax.set_xticks(major_ticks_x)
# ax.set_xticks(minor_ticks_x, minor=True)
# ax.set_yticks(major_ticks_y)
# ax.set_yticks(minor_ticks_y, minor=True)

markers = ['p','s','.','^','d','+','*',' ']
# markers = ['$a$','$b$','$c$','$d$','$e$','$f$','$g$','$h$']

xx_element        = np.array([1, 2, 3, 4, 5, 10, 15, 20])

klassisch_element = np.array([0.031749738, 0.032195229, 0.03216323, 0.03143856, 0.031160555, 0.031117369, 0.031145601, 0.031159884])
ad_element        = np.array([0.002113984, 0.002137648, 0.001997898, 0.001774646, 0.001647325, 0.001304078, 0.001174623, 0.001135134])
dif_element       = np.array([0.029635754, 0.030057581, 0.030165332, 0.029663914, 0.02951323, 0.029813291, 0.029970977, 0.03002475])


xx_polynom        = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10])
klassisch_polynom = np.array([0.032183421, 0.03218677, 0.032175384, 0.032657074, 0.033514649, 0.034861804, 0.036864568, 0.03966297, 0.043274208])
ad_polynom        = np.array([0.002129049, 0.00191033, 0.001629352, 0.003292059, 0.004072005, 0.00502966, 0.005869516, 0.006244993, 0.007005337])
dif_polynom       = np.array([0.030050578, 0.030277215, 0.030497683, 0.029308278, 0.029372675, 0.029743802, 0.031076943, 0.033447956, 0.03626053])

fig = plt.figure('Berechnungdauer Variation Elementanzahl')
ax = fig.add_subplot(1,1,1)
plt.tight_layout()

plt.xlabel('Anzahl Elemente')
plt.ylabel('Berechnungsdauer [s]')

plt.plot(xx_element, klassisch_element, linewidth=2.0, linestyle='-' , color='#005293' , label='klassisch Variation')
plt.plot(xx_element, ad_element, linewidth=2.0, linestyle='-' , color='#A2AD00' , label='hyper-duale Variation')
plt.plot(xx_element, dif_element, linewidth=2.0, linestyle='-' , color='black' , label='Differenz')

plt.legend(loc=9, bbox_to_anchor=(0.5, 1.2),
           ncol=3, borderaxespad=0.)

xx_lable = ['E1', 'E2', 'E3', 'E4', 'E5', 'E10', 'E15', 'E20']
major_ticks_x = np.array([1, 2, 3, 4, 5, 10, 15, 20])
major_ticks_y = np.arange(0,0.04,0.01)

ax.set_xticks(major_ticks_x)
ax.set_yticks(major_ticks_y)
ax.set_xticklabels(xx_lable)
plt.grid()

# matplotlib2tikz.save('D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/12_runtime_benchmark/runtime_element.tex',
#            figureheight = '\\figureheight',
#            figurewidth = '\\figurewidth')

# tikz_save('D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/12_runtime_benchmark/runtime_element.tikz',
#            figureheight = '\\figureheight',
#            figurewidth = '\\figurewidth')



fig = plt.figure('Berechnungdauer Variation Polynomgrad')
ax = fig.add_subplot(1,1,1)
plt.tight_layout()


plt.xlabel('Polynomgrad $p$ der Basisfunktionen')
plt.ylabel('Berechnungsdauer [s]')

plt.plot(xx_polynom, klassisch_polynom, linewidth=2.0, linestyle='-' , color='#005293' , label='klassisch Variation')
plt.plot(xx_polynom, ad_polynom, linewidth=2.0, linestyle='-' , color='#A2AD00' , label='hyper-duale Variation')
plt.plot(xx_polynom, dif_polynom, linewidth=2.0, linestyle='-' , color='black' ,   label='Differenz')

plt.legend(loc=9, bbox_to_anchor=(0.5, 1.2),
           ncol=3, borderaxespad=0.)

xx_lable = ['p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10']
major_ticks_x = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10])
major_ticks_y = np.arange(0,0.05,0.01)
ax.set_xticks(major_ticks_x)
ax.set_yticks(major_ticks_y)
ax.set_xticklabels(xx_lable)
# ax.set_aspect(aspect=60)

plt.grid()

tikz_save('D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/12_runtime_benchmark/runtime_degree.tikz',
           figureheight = '\\figureheight',
           figurewidth = '\\figurewidth')

# plt.legend(loc='best',  prop={'size':10}, bbox_to_anchor=(1.6, 1))
# plt.grid(alpha=0.5, linestyle='dotted')
# ax.set_yticks(major_ticks_y)

# # And a corresponding grid
# ax.grid(which='both')
# ax.grid(which='major', linestyle='dotted', alpha=0.8)





# # # plt.suptitle(r"$B$-spline basis functions, $p=%d$" % (3))
# plt.savefig("D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/kragarm_benchmark.pgf")
# plt.savefig("D:/OneDrive/Uni/Masterarbeit/Template_Thesis/images/10_Benchmark/kragarm_benchmark.pdf")

plt.show()



# plt.figure('Polyn. Grad',figsize=(6.6, 3.5))
# plt.clf()

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