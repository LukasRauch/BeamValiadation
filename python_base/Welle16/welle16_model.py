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
# import matplotlib as mpl
# import matplotlib.path as mpath
# import matplotlib.patches as mpatches
# import matplotlib.pyplot as plt
from geomdl import BSpline
from geomdl import utilities
from geomdl import exchange
from geomdl import Multi

from geomdl.visualization import VisMPL

start_time = time.time()

print('Process ID: ', os.getpid())
print(' ')

model = an.Model.open(r'C:\_Masterarbeit\beispiele\Welle16\Welle16_model_import.iga')

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
element_properties.SetValue(CROSS_AREA          , 5.3236E-03)     # m²
element_properties.SetValue(YOUNG_MODULUS       , 210000)      # kN/m²
element_properties.SetValue(SHEAR_MODULUS       , 81000)     # kN/m²
element_properties.SetValue(MOMENT_OF_INERTIA_Y , 1.143E-05)  # m4
element_properties.SetValue(MOMENT_OF_INERTIA_Z , 1.143E-05)  # m4
element_properties.SetValue(MOMENT_OF_INERTIA_T , 6.735E-07) # m4
element_properties.SetValue(POISSON_RATIO       , 0)        # m4

kratos_curve = NodeCurveGeometry3D(Degree = curve_geometry.Degree, NumberOfNodes = curve_geometry.NbPoles)

# Erstellen der Elemente
nodes = []
node_indices = []

for i in range(curve_geometry.NbPoles): # Erzeugung der 4 Kontrollpunkte mit den entsprechenden Koordinaten
                                    # ID                          X,                        Y,                         Z
    nodes = model_part.CreateNewNode(i+1, curve_geometry.Pole(i)[0], curve_geometry.Pole(i)[1], curve_geometry.Pole(i)[2])
    nodes.SetValue(NURBS_CONTROL_POINT_WEIGHT, curve_geometry.Weight(i))
    node_indices.append(nodes.Id)
    kratos_curve.SetNode(Index = i, Value = nodes)

# Erzeugen der Kurve
for i in range(curve_geometry.NbKnots):
    value = curve_geometry.Knot(i)
    kratos_curve.SetKnot(Index = i, Value = value)

print(kratos_curve.Knots())

# Berechnen der Formfunktionen an den entsprechnenden Gausspunktne
integration_points = curve_item.IntegrationPoints()
shapes = an.CurveShapeEvaluator(Degree = curve_geometry.Degree, Order = 2)

# Preprozessor Definitionen
t0 = [1, 0, 0]                  # Manuelle Vorgabe des Tangentenvektors
n0 = [0, 0, 1]                  # Manuelle Vorgabe des Normalenvektors
phi = 0                         # manuelle Vorgabe der Rotation
phi_der = 0                     # manuelle Vorgabe der Rotation 1st Ableitung
knots = curve_geometry.Knots

for n, (t, weight) in enumerate(integration_points):    # 4 Integrationspunkte

    # Formfunktionen an Gausspunkt n
    n_0 = Vector(4)                                     # Pro Integrationspunkt 4 Formfunktionauswertungen
    n_der = Matrix(2,4)                                 # Pro Integrationspunkt 2 x 4 Ableitungen
    shapes.Compute(curve_geometry.Knots, t)

    # ids der knoten die einen einfluss an der stelle t haben
    node_indices = [kratos_curve.Node(j).Id for j in
        range(shapes.FirstNonzeroPole, shapes.LastNonzeroPole + 1)]

    n_0 = [0] * shapes.NbNonzeroPoles
    n_1 = [0] * shapes.NbNonzeroPoles
    n_2 = [0] * shapes.NbNonzeroPoles

    for i in range(shapes.NbNonzeroPoles):
        n_0[i] = shapes(0, i)
        n_der[0,i] = shapes(1, i)
        n_der[1,i] = shapes(2, i)

    # Tangentenvektor ausgewertet an Gausspunkt n
    point, tangent = kratos_curve.DerivativesAt(T=t, Order=1)  # Tangentenvektor am aktuellen Integrationspunkt auswerten

    # Generierung der Elemente pro Integrationspunkt
    element = model_part.CreateNewElement('IgaBeamElement', n+1, node_indices, element_properties)
    element.SetValue(INTEGRATION_WEIGHT, weight)
    element.SetValue(SHAPE_FUNCTION_VALUES, n_0)                # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, n_der)   # Typ Matrix
    element.SetValue(KNOTSPANN, knots)                          # Knotspannvector der gesamten Kurve mitgeben
    element.SetValue(T0, tangent)

    ### mauelle Vorgbe
    element.SetValue(N0, n0)
    element.SetValue(PHI, phi)
    element.SetValue(PHI_0_DER, phi_der)

# Freiheitsgrade einfügen
dof_node = 4
VariableUtils().AddDof(DISPLACEMENT_X, REACTION_X, model_part)
VariableUtils().AddDof(DISPLACEMENT_Y, REACTION_Y, model_part)
VariableUtils().AddDof(DISPLACEMENT_Z, REACTION_Z, model_part)
VariableUtils().AddDof(DISPLACEMENT_ROTATION, REACTION_ROTATION, model_part)

# Randbedingungen: Auflager
# Kontrollpunkt 1
model_part.GetNode(1).Fix(DISPLACEMENT_X)
model_part.GetNode(1).Fix(DISPLACEMENT_Y)
model_part.GetNode(1).Fix(DISPLACEMENT_Z)
model_part.GetNode(1).Fix(DISPLACEMENT_ROTATION)

# model_part.GetNode(2).Fix(DISPLACEMENT_X)
model_part.GetNode(2).Fix(DISPLACEMENT_Y)
model_part.GetNode(2).Fix(DISPLACEMENT_Z)

model_part.GetNode(3).Fix(DISPLACEMENT_Y)
model_part.GetNode(4).Fix(DISPLACEMENT_Y)
model_part.GetNode(5).Fix(DISPLACEMENT_Y)
model_part.GetNode(6).Fix(DISPLACEMENT_Y)
model_part.GetNode(7).Fix(DISPLACEMENT_Y)
model_part.GetNode(8).Fix(DISPLACEMENT_Y)
model_part.GetNode(9).Fix(DISPLACEMENT_Y)
model_part.GetNode(10).Fix(DISPLACEMENT_Y)
model_part.GetNode(11).Fix(DISPLACEMENT_Y)
model_part.GetNode(12).Fix(DISPLACEMENT_Y)
model_part.GetNode(13).Fix(DISPLACEMENT_Y)
model_part.GetNode(14).Fix(DISPLACEMENT_Y)
model_part.GetNode(15).Fix(DISPLACEMENT_Y)
model_part.GetNode(16).Fix(DISPLACEMENT_Y)

# Randbedingungen: Knotenlast
load_properties = model_part.GetProperties()[2] # propperty-ID = 2
#                             typ,                     Id,  Knoten                   , Eigenschaften
model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [model_part.GetNode(16).Id], load_properties)


# Löser konfigurieren
model_part.SetBufferSize(3)

# Verfahren
time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

linear_solver = new_linear_solver_factory.ConstructSolver(Parameters(
    # r'{"solver_type": "SkylineLUFactorizationSolver"}'))
    r'{"solver_type": "eigen_sparse_lu"}'))

# Abbruchkriterium
relative_tolerance = 1e-7
absolute_tolerance = 1e-7
conv_criteria = ResidualCriteria(relative_tolerance, absolute_tolerance)
conv_criteria.SetEchoLevel(2)

# Löser
maximum_iterations = 3 #!! Wenn der Löser nur eine Iteration durchführt erhälst du eine lineare Lösung > Iterationszahl erhöhen!
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
solver.SetEchoLevel(0)

num_pole = curve_geometry.NbPoles
num_load_steps = 20  

disp_X = []
disp_Y = []
disp_Z = []

disp_X = np.empty([num_load_steps, num_pole])
disp_Y = np.empty([num_load_steps, num_pole])
disp_Z = np.empty([num_load_steps, num_pole])

for i in range(num_load_steps):
    F = i * 0.1/num_load_steps
    # node_2.SetSolutionStepValue(POINT_LOAD_Y, 1000 * (i + 1) / 10)
    model_part.GetNode(16).SetSolutionStepValue(POINT_LOAD_Z, F)

    # aktuellen modellzustand kopieren
    model_part.CloneTimeStep(i+1)

    # aktuellen zustand lösen
    print("solver step: ", i)
    solver.Solve()

    for j in range(curve_geometry.NbPoles):
        disp_X[i,j] = (model_part.GetNode(j+1).X )
        disp_Y[i,j] = (model_part.GetNode(j+1).Y )
        disp_Z[i,j] = (model_part.GetNode(j+1).Z )

print("\n\nStep " + str(i) + " :: Knoten Z: " + str(model_part.GetNode(num_pole).Z))
print("Verschiebung in X: " + str(model_part.GetNode(num_pole).X - model_part.GetNode(num_pole).X0))
print("Verschiebung in Y: " + str(model_part.GetNode(num_pole).Y - model_part.GetNode(num_pole).Y0))
print("Verschiebung in Z: " + str(model_part.GetNode(num_pole).Z - model_part.GetNode(num_pole).Z0))

print("Prozes time:  %s seconds ---" % (time.time() - start_time))
print('Calculations done!')


# Fix file path
# os.chdir(os.path.dirname(os.path.realpath(__file__)))
multi_curve = Multi.MultiCurve()

# Set up the Curve
for n in range(num_load_steps):
    curve = BSpline.Curve()
    curve.degree = 3

    curve.ctrlpts =[(disp_X[n,0], -disp_Z[n,0]),
                    (disp_X[n,1], -disp_Z[n,1]),
                    (disp_X[n,2], -disp_Z[n,2]),
                    (disp_X[n,3], -disp_Z[n,3]),
                    (disp_X[n,4], -disp_Z[n,4]),
                    (disp_X[n,5], -disp_Z[n,5]),
                    (disp_X[n,6], -disp_Z[n,6]),
                    (disp_X[n,7], -disp_Z[n,7]),
                    (disp_X[n,8], -disp_Z[n,8]),
                    (disp_X[n,9], -disp_Z[n,9]),
                    (disp_X[n,10], -disp_Z[n,10]),
                    (disp_X[n,11], -disp_Z[n,11]),
                    (disp_X[n,12], -disp_Z[n,12]),
                    (disp_X[n,13], -disp_Z[n,13]),
                    (disp_X[n,14], -disp_Z[n,14]),
                    (disp_X[n,15], -disp_Z[n,15])]
                    
    curve.knotvector = utilities.generate_knot_vector(curve.degree, len(curve.ctrlpts))
    # curve.knotvector = (0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 13.0, 13.0, 13.0)

    multi_curve.add(curve)

# curve0 = BSpline.Curve()
# curve1 = BSpline.Curve()
# curve2 = BSpline.Curve()
# curve0.degree = 3
# curve1.degree = 3
# curve2.degree = 3

# curve0.ctrlpts =[(disp_X[0,0], -disp_Z[0,0]),
#                 (disp_X[0,1], -disp_Z[0,1]),
#                 (disp_X[0,2], -disp_Z[0,2]),
#                 (disp_X[0,3], -disp_Z[0,3]),
#                 (disp_X[0,4], -disp_Z[0,4]),
#                 (disp_X[0,5], -disp_Z[0,5]),
#                 (disp_X[0,6], -disp_Z[0,6]),
#                 (disp_X[0,7], -disp_Z[0,7]),]
# curve0.knotvector = (0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0)
# curve1.ctrlpts =[(disp_X[1,0], -disp_Z[1,0]),
#                 (disp_X[1,1], -disp_Z[1,1]),
#                 (disp_X[1,2], -disp_Z[1,2]),
#                 (disp_X[1,3], -disp_Z[1,3]),
#                 (disp_X[1,4], -disp_Z[1,4]),
#                 (disp_X[1,5], -disp_Z[1,5]),
#                 (disp_X[1,6], -disp_Z[1,6]),
#                 (disp_X[1,7], -disp_Z[1,7]),]
# curve1.knotvector = (0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0)
# curve2.ctrlpts =[(disp_X[2,0], -disp_Z[2,0]),
#                 (disp_X[2,1], -disp_Z[2,1]),
#                 (disp_X[2,2], -disp_Z[2,2]),
#                 (disp_X[2,3], -disp_Z[2,3]),
#                 (disp_X[2,4], -disp_Z[2,4]),
#                 (disp_X[2,5], -disp_Z[2,5]),
#                 (disp_X[2,6], -disp_Z[2,6]),
#                 (disp_X[2,7], -disp_Z[2,7]),]
# curve2.knotvector = (0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0)

# multi_curve_x = Multi.MultiCurve(curve0, curve1, curve2)

multi_curve.delta = 0.05
vis_config = VisMPL.VisConfig(legend=False)
multi_curve.vis = VisMPL.VisCurve2D(vis_config)
multi_curve.render()

print("Prozes time:  %s seconds ---" % (time.time() - start_time))
print("Exit!")


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
