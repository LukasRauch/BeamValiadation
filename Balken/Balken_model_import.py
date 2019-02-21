# Module Importieren
from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import new_linear_solver_factory
import os
import yaml
import json
import time
import ANurbs as an
import numpy as np
import matplotlib as mpl
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from geomdl import BSpline
from geomdl import utilities
from geomdl import exchange
from geomdl import operations
from geomdl import Multi
from geomdl.visualization import VisMPL




start_time = time.time()

print('Process ID: ', os.getpid())
print(' ')

model = an.Model.open(r'C:\_Masterarbeit\BeamValidation\Balken\Balken.iga')

curve_item = model.of_type('Curve3D')[0]
curve = curve_item.data
curve_geometry = curve_item.geometry().geometry

#Modell Erzeugen
model_part = ModelPart('Model')

# variablen definieren die jeder knoten speichern soll
model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
model_part.AddNodalSolutionStepVariable(DISPLACEMENT_ROTATION)
model_part.AddNodalSolutionStepVariable(REACTION)
model_part.AddNodalSolutionStepVariable(REACTION_ROTATION)
model_part.AddNodalSolutionStepVariable(POINT_LOAD)
# model_part.AddNodalSolutionStepVariable(LOAD_VECTOR_MOMENT)

# Querschnittswerte
a = 1       # Querschnittshöhe
b = 1       # Querschnittsbreite
A = a * b   # Querschnittsfläche
Iy = (a * b^3)/12   # Flächenträgheitsmoment Iy
Iz = (b * a^3)/12   # Flächenträgheitsmoment Iz

# elementeigenschaften definieren
element_properties = model_part.GetProperties()[1] # property-id = 1
element_properties.SetValue(CROSS_AREA          , 100)     # m²
element_properties.SetValue(YOUNG_MODULUS       , 1)      # kN/m²
element_properties.SetValue(SHEAR_MODULUS       , 0.5)     # kN/m²
element_properties.SetValue(MOMENT_OF_INERTIA_Y , 100)  # m4
element_properties.SetValue(MOMENT_OF_INERTIA_Z , 500)  # m4
element_properties.SetValue(MOMENT_OF_INERTIA_T , 100) # m4
element_properties.SetValue(POISSON_RATIO       , 0)        # m4
element_properties.SetValue(DENSITY             , 78.5)

kratos_curve = NodeCurveGeometry3D(Degree = curve_geometry.Degree, NumberOfNodes = curve_geometry.NbPoles)

# Erstellen der Elemente
# nodes = []
node_indices = []
element_count = 0

for i in range(curve_geometry.NbPoles): # Erzeugung der 4 Kontrollpunkte mit den entsprechenden Koordinaten
                                    # ID                          X,                        Y,                         Z
    node = model_part.CreateNewNode(i+1, curve_geometry.Pole(i)[0], curve_geometry.Pole(i)[1], curve_geometry.Pole(i)[2])
    node.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
    # node.SetValue(NURBS_CONTROL_POINT_WEIGHT, curve_geometry.Weight(i))
    node_indices.append(node.Id)
    kratos_curve.SetNode(Index = i, Value = node)

# Erzeugen der Kurve
for i in range(curve_geometry.NbKnots):
    value = curve_geometry.Knot(i)
    kratos_curve.SetKnot(Index = i, Value = value)

# print(kratos_curve.Knots())

# Berechnen der Formfunktionen an den entsprechenden Gauss-punkte
integration_points = curve_item.IntegrationPoints()
shapes = an.CurveShapeEvaluator(Degree = curve_geometry.Degree, Order = 3)

# Preprocessor Definitionen
n0 = [0, 0, -1]                  # Manuelle Vorgabe des Normalenvektors
phi = 0                         # manuelle Vorgabe der Rotation
phi_der = 0                     # manuelle Vorgabe der Rotation 1st Ableitung
knots = curve_geometry.Knots

for n, (t, weight) in enumerate(integration_points):    # 4 Integrationspunkte
    element_count += 1

    # Geometie
    A2 = Vector(3)
    A2_1 = Vector(3)
    A3 = Vector(3) 
    A3_1 = Vector(3) 

    # Formfunktionen an Gausspunkt n
    n_0 = Vector(4)                                     
    n_1 = Vector(4)                                     
    n_2 = Vector(4)                                     
    n_3 = Vector(4)                                     
    # n_der = Matrix(2,4)                                 
    shapes.Compute(curve_geometry.Knots, t)

    for i in range(shapes.NbNonzeroPoles):
        n_0[i] = shapes(0, i)
        n_1[i] = shapes(1, i)
        n_2[i] = shapes(2, i)
        n_3[i] = shapes(3, i)

    # Tangentenvektor ausgewertet an Gausspunkt n
    # Normierung des Tangentenvektors erfolgt Kratos-intern
    point, tangent = kratos_curve.DerivativesAt(T=t, Order=1)  # Tangentenvektor am aktuellen Integrationspunkt auswerten

    # Generierung der Elemente pro Integrationspunkt
    # element = model_part.CreateNewElement('IgaBeamADElement', n+1, node_indices, element_properties)
    # element = model_part.CreateNewElement('IgaBeamWeakDirichletCondition', n+1, node_indices, element_properties)
    element = model_part.CreateNewElement('IgaBeamElement', n+1, node_indices, element_properties)
    element.SetValue(INTEGRATION_WEIGHT                 , weight)  # *2
    element.SetValue(SHAPE_FUNCTION_VALUES              , n_0)     # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DER_1         , n_1)     # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DER_2         , n_2)     # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DER_3         , n_3)     # Typ Vektor
    element.SetValue(T0                                 , tangent)
    element.SetValue(T0_DER                                 , [0,0,0])
    ### manuelle Vorgabe
    element.SetValue(N0                                 , n0)
    element.SetValue(PHI                                , phi)
    element.SetValue(PHI_DER_1                          , phi_der)


# Randbedingungen: Knotenlast
load_properties = model_part.GetProperties()[2] # property-ID = 2
#                             typ,                     Id,  Knoten                   , Eigenschaften
model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [model_part.GetNode(4).Id], load_properties)


# # _________________________________________________________________________________________________________________
# # Definition: Bettung
position_t = 0


# Formfunktionen an Gausspunkt n
n_0 = Vector(4)                                     
n_1 = Vector(4)                                     
n_2 = Vector(4)                                     
n_3 = Vector(4)                                     
# n_der = Matrix(2,4)                                 
shapes.Compute(curve_geometry.Knots, position_t)

for i in range(shapes.NbNonzeroPoles):
    n_0[i] = shapes(0, i)
    n_1[i] = shapes(1, i)
    n_2[i] = shapes(2, i)
    n_3[i] = shapes(3, i)

# Tangentenvektor ausgewertet an Gausspunkt n
# Normierung des Tangentenvektors erfolgt Kratos-intern
point, tangent = kratos_curve.DerivativesAt(T=position_t, Order=1)  # Tangentenvektor am aktuellen Integrationspunkt auswerten

# Generierung der Elemente pro Integrationspunkt
# element = model_part.CreateNewElement('IgaBeamADElement', n+1, node_indices, element_properties)
element_dirichlet_condition = model_part.CreateNewElement('IgaBeamWeakDirichletCondition', element_count+1, node_indices, element_properties)
element_dirichlet_condition.SetValue(INTEGRATION_WEIGHT                 , 1)  # *2
element_dirichlet_condition.SetValue(SHAPE_FUNCTION_VALUES              , n_0)     # Typ Vektor
element_dirichlet_condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_1         , n_1)     # Typ Vektor
element_dirichlet_condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_2         , n_2)     # Typ Vektor
element_dirichlet_condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_3         , n_3)     # Typ Vektor
element_dirichlet_condition.SetValue(T0                                 , tangent)
element_dirichlet_condition.SetValue(T0_DER                             , [0,0,0])
### manuelle Vorgabe
element_dirichlet_condition.SetValue(N0                                 , n0)
element_dirichlet_condition.SetValue(PHI                                , phi)
element_dirichlet_condition.SetValue(PHI_DER_1                          , phi_der)
### Randbedingungen 
element_dirichlet_condition.SetValue(PENALTY_DISPLACEMENT               , 1e12)
element_dirichlet_condition.SetValue(PENALTY_ROTATION                   , 1e12)
element_dirichlet_condition.SetValue(PENALTY_TORSION                    , 1e12)
element_dirichlet_condition.SetValue(DIRICHLET_CONDITION_TYPE           , 123)    # 1 Displacement, 2 Torsion , 3 Rotation Winkel, 4 Steigung

# # #_________________________________________________________________________________________________________________

# Freiheitsgrade einfügen
dof_node = 4
VariableUtils().AddDof(DISPLACEMENT_X, REACTION_X, model_part)
VariableUtils().AddDof(DISPLACEMENT_Y, REACTION_Y, model_part)
VariableUtils().AddDof(DISPLACEMENT_Z, REACTION_Z, model_part)
VariableUtils().AddDof(DISPLACEMENT_ROTATION, REACTION_ROTATION, model_part)

# Randbedingungen: Auflager
# Kontrollpunkt 1
# model_part.GetNode(1).Fix(DISPLACEMENT_X)
# model_part.GetNode(1).Fix(DISPLACEMENT_Y)
# model_part.GetNode(1).Fix(DISPLACEMENT_Z)
# model_part.GetNode(1).Fix(DISPLACEMENT_ROTATION)

# model_part.GetNode(2).Fix(DISPLACEMENT_X)
# model_part.GetNode(2).Fix(DISPLACEMENT_Y)
# model_part.GetNode(2).Fix(DISPLACEMENT_Z)
# model_part.GetNode(2).Fix(DISPLACEMENT_ROTATION)

# model_part.GetNode(3).Fix(DISPLACEMENT_Y)
# model_part.GetNode(4).Fix(DISPLACEMENT_Y)

# model_part.GetNode(3).Fix(DISPLACEMENT_ROTATION)
# model_part.GetNode(4).Fix(DISPLACEMENT_ROTATION)

# Löser konfigurieren
model_part.SetBufferSize(1)

# Verfahren
time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

linear_solver = new_linear_solver_factory.ConstructSolver(Parameters(
    r'{"solver_type": "eigen_sparse_lu"}'))
    # r'{"solver_type": "SkylineLUFactorizationSolver"}'))

# Abbruchkriterium
relative_tolerance = 1e-4
absolute_tolerance = 1e-4
conv_criteria = ResidualCriteria(relative_tolerance, absolute_tolerance)
conv_criteria.SetEchoLevel(1)

# Löser
maximum_iterations = 1000
compute_reactions = True
reform_dofs_at_each_iteration = True
move_mesh_flag = True

solver = ResidualBasedNewtonRaphsonStrategy(
    model_part,
    time_scheme,
    linear_solver,
    conv_criteria,
    maximum_iterations,
    compute_reactions,
    reform_dofs_at_each_iteration,
    move_mesh_flag
)
solver.SetEchoLevel(1)

num_pole = curve_geometry.NbPoles
num_load_steps = 5

disp_X = []
disp_Y = []
disp_Z = []

disp_X = np.empty([num_load_steps, num_pole])
disp_Y = np.empty([num_load_steps, num_pole])
disp_Z = np.empty([num_load_steps, num_pole])

# PLOT _________________________________________________________________________________
multi_curve = Multi.MultiCurve()

for i in range(0, num_load_steps+1):
    F = -i * 10/num_load_steps
    moment_vec          = [0, i * 10/num_load_steps, 0]
    model_part.GetNode(4).SetSolutionStepValue(POINT_LOAD_Z, F)
    # model_part.GetElement(5)
    # element_load_properties.SetValue(LOAD_VECTOR_MOMENT, moment_vec)


    # aktuellen modellzustand kopieren
    model_part.CloneTimeStep(i+1)

    # aktuellen zustand lösen
    #____________________________________________________________________________________
    #####################################################################################
    #####################################################################################
    print("\nsolver step: ", i, "F =", F)
    solver.Solve()
    #____________________________________________________________________________________



    print("\n  COORDINATES ")
    print("DOF-TYPE" + "\t" + "X" + "\t\t\t" +  "Y" + "\t\t\t" + "Z")
    print("# Nr. =============================================================================")

    for k in range(curve_geometry.NbPoles):
        print("Pole "+ str(k+1) + "\t\t"
                    + '%.12f' % (model_part.GetNode(k+1).X ) +  "\t\t"
                    + '%.12f' % (model_part.GetNode(k+1).Y ) +  "\t\t"
                    + '%.12f' % (model_part.GetNode(k+1).Z ) +  "\t\t" )

    print("\n  DISPLACEMENTS ")
    print("DOF-TYPE" + "\t" + "X" + "\t\t\t" +  "Y" + "\t\t\t" + "Z")
    print("# Nr. =============================================================================")

    for k in range(curve_geometry.NbPoles):
        print("Pole "+ str(k+1) + "\t\t"
                    + '%.12f' % (model_part.GetNode(k+1).X - model_part.GetNode(k+1).X0) +  "\t\t"
                    + '%.12f' % (model_part.GetNode(k+1).Y - model_part.GetNode(k+1).Y0) +  "\t\t"
                    + '%.12f' % (model_part.GetNode(k+1).Z - model_part.GetNode(k+1).Z0) +  "\t\t" )


    # Create a B-Spline Curve instance
    plot_curve = BSpline.Curve()
    #Set up a Curve
    plot_curve.degree = curve_geometry.Degree
    ctlpts_list = open("C:\_Masterarbeit\BeamValidation\Balken\ctlpts_list.txt", "w")

    # Draw the control points polygon, the 3D curve and the vectors
    # fig = plt.figure('figure X', figsize=(10.67, 8), dpi= 96) 
    # ax = Axes3D(fig)
    
    for j in range(curve_geometry.NbPoles):
        data = model_part.GetNode(j+1).X , model_part.GetNode(j+1).Y , model_part.GetNode(j+1).Z
        ctlpts_list.write(str(model_part.GetNode(j+1).X ) + ',' +
                          str(model_part.GetNode(j+1).Y ) + ',' +
                          str(model_part.GetNode(j+1).Z ) + '\n')

    ctlpts_list.close()
    plot_curve.ctrlpts = exchange.import_txt("C:\_Masterarbeit\BeamValidation\Balken\ctlpts_list.txt")
    plot_curve.knotvector = utilities.generate_knot_vector(plot_curve.degree, len(plot_curve.ctrlpts))
    # Set evaluation delta 
    # plot_curve.delta = 0.001
    # plot_curve.evaluate 

    # add to mulit_curve
    multi_curve.add(plot_curve)


print("\nProzes time:  %s seconds ---" % (time.time() - start_time))
print('done!')
# print("Prozes time:  %s seconds ---" % (time.time() - start_time))
# print("Exit!")print("\nProzes time:  %s seconds ---" % (time.time() - start_time))
print('done!')
# print("Prozes time:  %s seconds ---" % (time.time() - start_time))
# print("Exit!")

# Set evaluation delta 
multi_curve.delta = 0.001
# multi_curve.evaluate
# plot the controlpoint polygon and the evaluated curve
vis_comp = VisMPL.VisCurve3D()
multi_curve.vis = vis_comp
multi_curve.render(cpcolor='black', evalcolor='red') 
pass




    
    
    
#     # List of parametric coordinates to be evaluated
#     u_list = curve_geometry.Knots
#     # Normalize u_list
#     u_list[:] = [x / u_list[-1] for x in u_list]

#     # Evaluate tangents, normals and binormals, respectively
#     curvetans       = [[] for _ in range(len(u_list))]
#     curvenorms      = [[] for _ in range(len(u_list))]
#     curvebionorms   = [[] for _ in range(len(u_list))]
#     for idx, u in enumerate(u_list):
#         curvetans[idx]      = operations.tangent(plot_curve, u, normalize=True)
#         curvenorms[idx]     = operations.normal(plot_curve, u, normalize=True)
#         curvebionorms[idx]  = operations.binormal(plot_curve, u, normalize=True)

    
#     # Control Points, Curve and Tangent Vector Plotting using Matplotlib 
#     ###

#     # Arrange control points and evaluate curve points for plotting
#     ctrlpts = np.array(plot_curve.ctrlpts)
#     curvepts = np.array(plot_curve.evalpts)

#     # convert tangent, normal and binormal vector lists into NumPy array
#     ctarr = np.array(curvetans)
#     cnarr = np.array(curvenorms)
#     cbnarr = np.array(curvebionorms)

#     # Plot the 3D lines
#     ax.plot(ctrlpts[:,0], ctrlpts[:,1], ctrlpts[:,2], color='black', linestyle='-', marker='o', linewidth=1)
#     ax.plot(curvepts[:, 0], curvepts[:, 1], curvepts[:, 2], color='red', linestyle='-', linewidth=2)

#     # # Plot tangent vectors
#     # ax.quiver(ctarr[:, 0, 0], ctarr[:, 0, 1], ctarr[:, 0, 2], ctarr[:, 1, 0], ctarr[:, 1, 1], ctarr[:, 1, 2],
#     #         color='blue', length=2.5)

#     # # Plot normal vectors
#     # ax.quiver(cnarr[:, 0, 0], cnarr[:, 0, 1], cnarr[:, 0, 2], cnarr[:, 1, 0], cnarr[:, 1, 1], cnarr[:, 1, 2],
#     #         color='red', length=2.5)

#     # # Plot binormal vectors
#     # ax.quiver(cbnarr[:, 0, 0], cbnarr[:, 0, 1], cbnarr[:, 0, 2], cbnarr[:, 1, 0], cbnarr[:, 1, 1], cbnarr[:, 1, 2],
#     #         color='green', length=2.5)

#     # # Add legend to 3D plot, @ref: https://stackoverflow.com/a/20505720
#     # ctrlpts_proxy = mpl.lines.Line2D([0], [0], linestyle='-.', color='black', marker='o')
#     # curvepts_proxy = mpl.lines.Line2D([0], [0], linestyle='none', color='brown', marker='o')
#     # tangent_proxy = mpl.lines.Line2D([0], [0], linestyle='none', color='blue', marker='>')
#     # normal_proxy = mpl.lines.Line2D([0], [0], linestyle='none', color='red', marker='>')
#     # binormal_proxy = mpl.lines.Line2D([0], [0], linestyle='none', color='green', marker='>')
#     # ax.legend([ctrlpts_proxy, curvepts_proxy, tangent_proxy, normal_proxy, binormal_proxy],
#     #         ['Control Points', 'Curve', 'Tangents', 'Normals', 'Binormals'], numpoints=1)

# # Display the 3D plot
# plt.show()



# Path = mpath.Path
# fig, ax = plt.subplots()

# # control_points = np.empty()

# for n in range(num_load_steps):
#     control_points =[(disp_X[n,0], -disp_Z[n,0]),
#                     (disp_X[n,1], -disp_Z[n,1]),
#                     (disp_X[n,2], -disp_Z[n,2]),
#                     (disp_X[n,3], -disp_Z[n,3])]

#     crv = mpatches.PathPatch(Path(control_points,
#     [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
#     fc = 'none',
#     transform = ax.transData)

#     ax.add_patch(crv)
#     ax.plot(disp_X[n,:], -disp_Z[n,:], 'rx',disp_X[n,:], -disp_Z[n,:], 'c:' )

# plt.show()
# plt.grid(True)
# plt.axis([-3,12,-20,0.5])







# print('Kontrollpunkte und Gewichte:')

# for i in range(curve_geometry.NbPoles):
#     pole = curve_geometry.Pole(i)
#     weight = curve_geometry.Weight(i)

#     print('  ', i, pole, weight)

# print()
# print('Knots:')
# knots = curve_geometry.Knots
# print('  ', knots)
# print()
# print('Integrationspunkte und Gewichte im Parameterraum der Kurve:')

# integration_points = curve_item.IntegrationPoints()

# for t, weight in integration_points:
#     print('  ', 't =', t, 'weight =', weight)

# print()
# print('Tangentenvektor:')

# integration_points = curve_item.IntegrationPoints()

# for t, weight in integration_points:
#     point, tangent = curve.DerivativesAt(T=t, Order=1)
#     print('  ', 't =', t, 'tangent =', tangent)

# print()
# print('Formfunktionen:')

# shapes = an.CurveShapeEvaluator(Degree=curve_geometry.Degree, Order=2)

# for t, weight in integration_points:
#     shapes.Compute(curve_geometry.Knots, t)

#     n_0 = [0] * shapes.NbNonzeroPoles
#     n_1 = [0] * shapes.NbNonzeroPoles
#     n_2 = [0] * shapes.NbNonzeroPoles

#     for i in range(shapes.NbNonzeroPoles):
#         n_0[i] = shapes(0, i)
#         n_1[i] = shapes(1, i)
#         n_2[i] = shapes(2, i)

#     print('  ', 't =', t)
#     print('    ', 'n_0 =', n_0)
#     print('    ', 'n_1 =', n_1)
#     print('    ', 'n_2 =', n_2)

# print(curve_geometry.PointAt(1.0))


#_________________________________________________________________________________________________________________
# # Definition: Moment
# moment_vec          = [0, 1, 0]
# moment_pos_para     = 10.0

# # # Formfunktionen an Gausspunkt n
# # n_0 = Vector(4)                                     
# # n_der = Matrix(2,4)                                 
# # shapes.Compute(curve_geometry.Knots, moment_pos_para)

# # for i in range(shapes.NbNonzeroPoles):
# #     n_0[i] = shapes(0, i)
# #     n_der[0,i] = shapes(1, i)
# #     n_der[1,i] = shapes(2, i)

# # Formfunktionen an Gausspunkt n
# n_0 = Vector(4)                                     
# n_1 = Vector(4)                                     
# n_2 = Vector(4)                                     
# n_3 = Vector(4)                                     
# # n_der = Matrix(2,4)                                 
# shapes.Compute(curve_geometry.Knots, moment_pos_para)

# for i in range(shapes.NbNonzeroPoles):
#     n_0[i] = shapes(0, i)
#     n_1[i] = shapes(1, i)
#     n_2[i] = shapes(2, i)
#     n_3[i] = shapes(3, i)

# # Tangentenvektor ausgewertet an Gausspunkt n
# # Normierung des Tangentenvektors erfolgt Kratos-intern
# point, tangent = kratos_curve.DerivativesAt(T=t, Order=1)  # Tangentenvektor am aktuellen Integrationspunkt auswerten

# element_load_properties = model_part.GetProperties()[3] # property-id = 3
# # element_load_properties.SetValue(LOAD_VECTOR_MOMENT                 , moment_vec)
# # Generierung der Elemente pro Integrationspunkt
# load_condition_element = model_part.CreateNewElement('IgaBeamMomentCondition', 5, node_indices, element_load_properties)
# # element_load_properties.SetValue(LOAD_VECTOR_MOMENT                  , moment_vec)
# load_condition_element.SetValue(INTEGRATION_WEIGHT                  , 1)  # *2
# load_condition_element.SetValue(SHAPE_FUNCTION_VALUES              , n_0)     # Typ Vektor
# load_condition_element.SetValue(SHAPE_FUNCTION_LOCAL_DER_1         , n_1)     # Typ Vektor
# load_condition_element.SetValue(SHAPE_FUNCTION_LOCAL_DER_2         , n_2)     # Typ Vektor
# load_condition_element.SetValue(SHAPE_FUNCTION_LOCAL_DER_3         , n_3)     # Typ Vektor
# load_condition_element.SetValue(T0                                 , tangent)

# ### manuelle Vorgabe
# load_condition_element.SetValue(N0, n0)
# load_condition_element.SetValue(PHI, phi)
# load_condition_element.SetValue(PHI_DER_1, phi_der)