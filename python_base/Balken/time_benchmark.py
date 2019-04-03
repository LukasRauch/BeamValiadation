# Module Importieren
from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import new_linear_solver_factory
import os
import ANurbs as an
import numpy as np
import time


start_time = time.time()

print('Process ID: ', os.getpid())
print(' ')

# model = an.Model.open(r'C:\_Masterarbeit\BeamValidation\Sweep\sweep.iga')
model = an.Model.open(r'C:\_Masterarbeit\BeamValidation\pathon_anurbs\data\Benchmark_Time\D3E5.iga')

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
element_properties.SetValue(CROSS_AREA          , 200)     # m²
element_properties.SetValue(YOUNG_MODULUS       , 250)      # kN/m²
element_properties.SetValue(SHEAR_MODULUS       , 175)     # kN/m²
element_properties.SetValue(MOMENT_OF_INERTIA_Y , 100)  # m4
element_properties.SetValue(MOMENT_OF_INERTIA_Z , 100)  # m4
element_properties.SetValue(MOMENT_OF_INERTIA_T , 50) # m4

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
n0 = [0, 0, 1]                  # Manuelle Vorgabe des Normalenvektors
t0_1 = [0, 0, 0]
phi = 0                         # manuelle Vorgabe der Rotation
phi_der = 0                     # manuelle Vorgabe der Rotation 1st Ableitung
knots = curve_geometry.Knots

for n, (t, weight) in enumerate(integration_points):    # 11 Integrationspunkte
    # Geometie
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
    element.SetValue(T0_DER                             , t0_1)

   ### manuelle Vor
    element.SetValue(N0                                 , n0)
    element.SetValue(PHI                                , phi)
    element.SetValue(PHI_DER_1                          , phi_der)

# Randbedingungen: Knotenlast
load_properties = model_part.GetProperties()[2] # propperty-ID = 2
#                             typ,                     Id,  Knoten                   , Eigenschaften
model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [model_part.GetNode(curve_geometry.NbPoles).Id], load_properties)

# # Freiheitsgrade einfügen
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
# Kontrollpunkt 2
model_part.GetNode(2).Fix(DISPLACEMENT_Y)
model_part.GetNode(2).Fix(DISPLACEMENT_Z)

# Löser konfigurieren
model_part.SetBufferSize(1)

# Verfahren
time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

linear_solver = new_linear_solver_factory.ConstructSolver(Parameters(
    # r'{"solver_type": "SkylineLUFactorizationSolver"}'))
    r'{"solver_type": "eigen_sparse_lu"}'))

# Abbruchkriterium
relative_tolerance = 1e-4
absolute_tolerance = 1e-4
conv_criteria = ResidualCriteria(relative_tolerance, absolute_tolerance)
conv_criteria.SetEchoLevel(1)

# Löser
maximum_iterations = 100 #!! Wenn der Löser nur eine Iteration durchführt erhälst du eine lineare Lösung > Iterationszahl erhöhen!
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

num_load_steps = 10

for i in range(1, num_load_steps+1):
    F           = i * 100/num_load_steps
    model_part.GetNode(curve_geometry.NbPoles).SetSolutionStepValue(POINT_LOAD_Z, F)
    
    # aktuellen modellzustand kopieren
    model_part.CloneTimeStep(i+1)

    # aktuellen zustand lösen
    print("Kraft F: ", F)
    print("solver step: ", i)
    solver.Solve()

print("done! ")

print("\n  DISPLACEMENTS ")
print("DOF-TYPE" + "\t" + "X" + "\t\t\t" +  "Y" + "\t\t\t" + "Z")
print("# Nr. =============================================================================")

for k in range(curve_geometry.NbPoles):
    print("Pole "+ str(k+1) + "\t\t"
                + '%.12f' % (model_part.GetNode(k+1).X - model_part.GetNode(k+1).X0) +  "\t\t"
                + '%.12f' % (model_part.GetNode(k+1).Y - model_part.GetNode(k+1).Y0) +  "\t\t"
                + '%.12f' % (model_part.GetNode(k+1).Z - model_part.GetNode(k+1).Z0) +  "\t\t" )