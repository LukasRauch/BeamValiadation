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

model = an.Model.open(r'C:\_Masterarbeit\BeamValidation\Balken8\Balken8.iga')

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

# elementeigenschaften definieren
element_properties = model_part.GetProperties()[1] # property-id = 1
element_properties.SetValue(CROSS_AREA          , 1000)      # m²
element_properties.SetValue(YOUNG_MODULUS       , 1)      # kN/m²
element_properties.SetValue(SHEAR_MODULUS       , 0.5)      # kN/m²
element_properties.SetValue(MOMENT_OF_INERTIA_Y , 100)      # m4
element_properties.SetValue(MOMENT_OF_INERTIA_Z , 500)      # m4
element_properties.SetValue(MOMENT_OF_INERTIA_T , 100)      # m4
element_properties.SetValue(POISSON_RATIO       , 0)                # m4
element_properties.SetValue(DENSITY             , 78.5)
kratos_curve = NodeCurveGeometry3D(Degree = curve_geometry.Degree, NumberOfNodes = curve_geometry.NbPoles)

# Erstellen der Elemente
node = []
node_indices = []

for i in range(curve_geometry.NbPoles): # Erzeugung der 4 Kontrollpunkte mit den entsprechenden Koordinaten
                                    # ID                          X,                        Y,                         Z
    node = model_part.CreateNewNode(i+1, curve_geometry.Pole(i)[0], curve_geometry.Pole(i)[1], curve_geometry.Pole(i)[2])
    node.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
    node_indices.append(node.Id)
    kratos_curve.SetNode(Index = i, Value = node)

# Erzeugen der Kurve
for i in range(curve_geometry.NbKnots):
    value = curve_geometry.Knot(i)
    kratos_curve.SetKnot(Index = i, Value = value)

# Berechnen der Formfunktionen an den entsprechnenden Gausspunktne
integration_points = curve_item.IntegrationPoints()
shapes = an.CurveShapeEvaluator(Degree = curve_geometry.Degree, Order = 3)

# Preprozessor Definitionen
n0 = [0, 0, -1]                  # Manuelle Vorgabe des Normalenvektors
phi = 0                         # manuelle Vorgabe der Rotation
phi_der = 0                     # manuelle Vorgabe der Rotation 1st Ableitung
knots = curve_geometry.Knots

element_count = 0

for n, (t, weight) in enumerate(integration_points):    # 11 Integrationspunkte
    element_count += 1
    n_0 = [0] * shapes.NbNonzeroPoles
    n_1 = [0] * shapes.NbNonzeroPoles
    n_2 = [0] * shapes.NbNonzeroPoles
    n_3 = [0] * shapes.NbNonzeroPoles

    # Formfunktionen an Gausspunkt n
    shapes.Compute(curve_geometry.Knots, t)

    # ids der knoten die einen einfluss an der stelle t haben
    node_indices = [kratos_curve.Node(j).Id for j in
        range(shapes.FirstNonzeroPole, shapes.LastNonzeroPole + 1)]

    # node_indices = kratos_curve.Node(j).Id for j in 

    for i in range(shapes.NbNonzeroPoles):
        n_0[i] = shapes(0, i)
        n_1[i] = shapes(1, i)
        n_2[i] = shapes(2, i)
        n_3[i] = shapes(3, i)

    # Tangentenvektor ausgewertet an Gausspunkt n
    point, tangent = kratos_curve.DerivativesAt(T=t, Order=1)  # Tangentenvektor am aktuellen Integrationspunkt auswerten

    # Generierung der Elemente pro Integrationspunkt
    element = model_part.CreateNewElement('IgaBeamElement', n+1, node_indices, element_properties)
    element.SetValue(INTEGRATION_WEIGHT, weight)
    element.SetValue(SHAPE_FUNCTION_VALUES              , n_0)     # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DER_1         , n_1)     # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DER_2         , n_2)     # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DER_3         , n_3)     # Typ Vektor
    element.SetValue(T0                                 , tangent)
    ### manuelle Vorgabe
    element.SetValue(N0                                 , n0)
    element.SetValue(PHI                                , phi)
    element.SetValue(PHI_DER_1                         , phi_der)

    # Randbedingungen: Knotenlast
load_properties = model_part.GetProperties()[2] # propperty-ID = 2
#                             typ,                     Id,  Knoten                   , Eigenschaften
model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [model_part.GetNode(8).Id], load_properties)

# # _________________________________________________________________________________________________________________
# # Definition: Bettung
position_t = 0


# Formfunktionen an Gausspunkt n
n_0 = [0] * shapes.NbNonzeroPoles
n_1 = [0] * shapes.NbNonzeroPoles
n_2 = [0] * shapes.NbNonzeroPoles
n_3 = [0] * shapes.NbNonzeroPoles                                    
# n_der = Matrix(2,4)                                 
shapes.Compute(curve_geometry.Knots, position_t)
node_indices = np.arange(shapes.FirstNonzeroPole+1, shapes.LastNonzeroPole+2, dtype=int)

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
# # 
# # Freiheitsgrade einfügen
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

# # model_part.GetNode(2).Fix(DISPLACEMENT_X)
# model_part.GetNode(2).Fix(DISPLACEMENT_Y)
# model_part.GetNode(2).Fix(DISPLACEMENT_Z)

# model_part.GetNode(3).Fix(DISPLACEMENT_Y)
# model_part.GetNode(4).Fix(DISPLACEMENT_Y)
# model_part.GetNode(5).Fix(DISPLACEMENT_Y)
# model_part.GetNode(6).Fix(DISPLACEMENT_Y)
# model_part.GetNode(7).Fix(DISPLACEMENT_Y)
# model_part.GetNode(8).Fix(DISPLACEMENT_Y)

# Löser konfigurieren
model_part.SetBufferSize(1)

# Verfahren
time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

linear_solver = new_linear_solver_factory.ConstructSolver(Parameters(
    # r'{"solver_type": "SkylineLUFactorizationSolver"}'))
    r'{"solver_type": "eigen_sparse_lu"}'))

# Abbruchkriterium
relative_tolerance = 1e-7
absolute_tolerance = 1e-7
conv_criteria = ResidualCriteria(relative_tolerance, absolute_tolerance)
conv_criteria.SetEchoLevel(1)

# Löser
maximum_iterations = 200 #!! Wenn der Löser nur eine Iteration durchführt erhälst du eine lineare Lösung > Iterationszahl erhöhen!
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
    model_part.GetNode(8).SetSolutionStepValue(POINT_LOAD_Z, F)
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