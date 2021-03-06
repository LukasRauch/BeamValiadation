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

    # Formfunktionen an Gausspunkt n
    n_0 = Vector(4)                                     # Pro Integrationspunkt 4 Formfunktionauswertungen
    n_1 = Vector(4)                                     # Pro Integrationspunkt 4 Formfunktionauswertungen
    n_2 = Vector(4)                                     # Pro Integrationspunkt 4 Formfunktionauswertungen
    n_3 = Vector(4)                                     # Pro Integrationspunkt 4 Formfunktionauswertungen
    # n_der = Matrix(2,4)                                 # Pro Integrationspunkt 2 x 4 Ableitungen
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
    element = model_part.CreateNewElement('IgaBeamElement', n+1, node_indices, element_properties)
    element.SetValue(INTEGRATION_WEIGHT                 , weight)  # *2
    element.SetValue(SHAPE_FUNCTION_VALUES              , n_0)     # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DER_1         , n_1)     # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DER_2         , n_2)     # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DER_3         , n_3)     # Typ Vektor
    element.SetValue(T0                                 , tangent)
    ### manuelle Vorgabe
    element.SetValue(N0                                 , n0)
    element.SetValue(PHI                                , phi)
    element.SetValue(PHI_0_DER                          , phi_der)


# Randbedingungen: Knotenlast
load_properties = model_part.GetProperties()[2] # property-ID = 2
#                             typ,                     Id,  Knoten                   , Eigenschaften
model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [model_part.GetNode(8).Id], load_properties)

#_________________________________________________________________________________________________________________
# # Definition: Moment 
# moment_vec          = [0, 1, 0]
# moment_pos_para     = 10

# # # Formfunktionen an Gausspunkt n
# # n_0 = Vector(4)                                     # Pro Integrationspunkt 4 Formfunktionauswertungen
# # n_der = Matrix(2,4)                                 # Pro Integrationspunkt 2 x 4 Ableitungen
# # shapes.Compute(curve_geometry.Knots, moment_pos_para)

# # for i in range(shapes.NbNonzeroPoles):
# #     n_0[i] = shapes(0, i)
# #     n_der[0,i] = shapes(1, i)
# #     n_der[1,i] = shapes(2, i)

# # Formfunktionen an Gausspunkt n
# n_0 = Vector(4)                                     # Pro Integrationspunkt 4 Formfunktionauswertungen
# n_1 = Vector(4)                                     # Pro Integrationspunkt 4 Formfunktionauswertungen
# n_2 = Vector(4)                                     # Pro Integrationspunkt 4 Formfunktionauswertungen
# n_3 = Vector(4)                                     # Pro Integrationspunkt 4 Formfunktionauswertungen
# # n_der = Matrix(2,4)                                 # Pro Integrationspunkt 2 x 4 Ableitungen
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
# load_condition_element = model_part.CreateNewElement('IgaBeamMomentCondition', n+2, node_indices, element_load_properties)
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
# load_condition_element.SetValue(PHI_0_DER, phi_der)  
#_________________________________________________________________________________________________________________

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
# model_part.GetNode(2).Fix(DISPLACEMENT_X)
model_part.GetNode(2).Fix(DISPLACEMENT_Y)
model_part.GetNode(2).Fix(DISPLACEMENT_Z)

model_part.GetNode(3).Fix(DISPLACEMENT_Y)
model_part.GetNode(4).Fix(DISPLACEMENT_Y)

# Löser konfigurieren
model_part.SetBufferSize(1)

# Verfahren
time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

linear_solver = new_linear_solver_factory.ConstructSolver(Parameters(
    r'{"solver_type": "eigen_sparse_lu"}'))
    # r'{"solver_type": "SkylineLUFactorizationSolver"}'))

# Abbruchkriterium
relative_tolerance = 1e-7
absolute_tolerance = 1e-7
conv_criteria = ResidualCriteria(relative_tolerance, absolute_tolerance)
conv_criteria.SetEchoLevel(1)

# Löser
maximum_iterations = 200
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
num_load_steps = 10

disp_X = []
disp_Y = []
disp_Z = []

disp_X = np.empty([num_load_steps, num_pole])
disp_Y = np.empty([num_load_steps, num_pole])
disp_Z = np.empty([num_load_steps, num_pole])

for i in range(1, num_load_steps+1):
    F = i * 1/num_load_steps
    moment_vec          = [0, i * 2/num_load_steps, 0]
    model_part.GetNode(8).SetSolutionStepValue(POINT_LOAD_Z, F)
    # model_part.GetElement(5)
    # element_load_properties.SetValue(LOAD_VECTOR_MOMENT, moment_vec)

    # aktuellen modellzustand kopieren
    model_part.CloneTimeStep(i+1)

    # aktuellen zustand lösen
    # print("solver step: ", i)
    print("\nsolver step: ", i, "F =", F)
    # print("Verhältnis F*L^2/EI= ", F/(element_properties.GetValue(YOUNG_MODULUS)*element_properties.GetValue(MOMENT_OF_INERTIA_Y) ))
    solver.Solve()

    for j in range(curve_geometry.NbPoles):
        disp_X[i-1,j] = (model_part.GetNode(j+1).X )
        disp_Y[i-1,j] = (model_part.GetNode(j+1).Y )
        disp_Z[i-1,j] = (model_part.GetNode(j+1).Z )
    


    print("\nDOF-TYPE" + "\t" + "X" + "\t\t\t" +  "Y" + "\t\t\t" + "Z")
    print("# Nr. =============================================================================")

    for k in range(curve_geometry.NbPoles):
        print("Pole "+ str(k+1) + "\t\t"  
                    + '%.12f' % (model_part.GetNode(k+1).X - model_part.GetNode(k+1).X0) +  "\t\t"
                    + '%.12f' % (model_part.GetNode(k+1).Y - model_part.GetNode(k+1).Y0) +  "\t\t" 
                    + '%.12f' % (model_part.GetNode(k+1).Z - model_part.GetNode(k+1).Z0) +  "\t\t" )
    


print("\nProzes time:  %s seconds ---" % (time.time() - start_time))
print('done!')

Path = mpath.Path
fig, ax = plt.subplots()

# control_points = np.empty()

for n in range(num_load_steps):
    control_points =[(disp_X[n,0], -disp_Z[n,0]),
                    (disp_X[n,1], -disp_Z[n,1]),
                    (disp_X[n,2], -disp_Z[n,2]),
                    (disp_X[n,3], -disp_Z[n,3])]

    crv = mpatches.PathPatch(Path(control_points, 
    [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]), 
    fc = 'none',
    transform = ax.transData)

    ax.add_patch(crv)
    ax.plot(disp_X[n,:], -disp_Z[n,:], 'rx',disp_X[n,:], -disp_Z[n,:], 'c:' )

plt.show()
plt.grid(True)
plt.axis([-3,12,-20,0.5])



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
