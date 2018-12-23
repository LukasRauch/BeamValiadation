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

start_time = time.time()

print('Process ID: ', os.getpid())


model = an.Model.open(r'C:\_Masterarbeit\beispiele\Balken\Balken.iga')

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
element_properties.SetValue(CROSS_AREA          , 0.09)
element_properties.SetValue(YOUNG_MODULUS       , 210000000)
element_properties.SetValue(SHEAR_MODULUS       , 79300000)
element_properties.SetValue(MOMENT_OF_INERTIA_Y , 0.000675)
element_properties.SetValue(MOMENT_OF_INERTIA_Z , 0.000675)
element_properties.SetValue(MOMENT_OF_INERTIA_T , 0.001134)

# Erstellen der Elemente
nodes = []
node_indices = []

node_1 = model_part.CreateNewNode(1,curve_geometry.Pole(0)[0], curve_geometry.Pole(0)[1], curve_geometry.Pole(0)[2] )
node_2 = model_part.CreateNewNode(2,curve_geometry.Pole(1)[0], curve_geometry.Pole(1)[1], curve_geometry.Pole(1)[2] )
node_3 = model_part.CreateNewNode(3,curve_geometry.Pole(2)[0], curve_geometry.Pole(2)[1], curve_geometry.Pole(2)[2] )
node_4 = model_part.CreateNewNode(4,curve_geometry.Pole(3)[0], curve_geometry.Pole(3)[1], curve_geometry.Pole(3)[2] )

node_1.SetValue(NURBS_CONTROL_POINT_WEIGHT, curve_geometry.Weight(0))
node_1.SetValue(NURBS_CONTROL_POINT_WEIGHT, curve_geometry.Weight(1))
node_1.SetValue(NURBS_CONTROL_POINT_WEIGHT, curve_geometry.Weight(2))
node_1.SetValue(NURBS_CONTROL_POINT_WEIGHT, curve_geometry.Weight(3))

curve = NodeCurveGeometry3D(Degree = curve_geometry.Degree, NumberOfNodes = curve_geometry.NbPoles)

curve.SetKnot(Index = 0, Value = 0)
curve.SetKnot(Index = 1, Value = 1)
curve.SetKnot(Index = 2, Value = 2)
curve.SetKnot(Index = 3, Value = 3)

curve.SetNode(Index = 0, Value = node_1)
curve.SetNode(Index = 1, Value = node_2)
curve.SetNode(Index = 2, Value = node_3)
curve.SetNode(Index = 3, Value = node_4)


# for i in range(curve_geometry.NbPoles): # Erzeugung der 4 Kontrollpunkte mit den entsprechenden Koordinaten
                                    # ID                          X,                        Y,                         Z
    # nodes = model_part.CreateNewNode(i+1, curve_geometry.Pole(i)[0], curve_geometry.Pole(i)[1], curve_geometry.Pole(i)[2])
    # nodes.SetValue(NURBS_CONTROL_POINT_WEIGHT, curve_geometry.Weight(i))  
    # node_indices.append(nodes.Id)

# Berechnen der Formfunktionen an den entsprechnenden Gausspunktne 
integration_points = curve_item.IntegrationPoints()
shapes = an.CurveShapeEvaluator(Degree = curve_geometry.Degree, Order = 2)


for n, (t, weight) in enumerate(integration_points):
    
    n_0 = Vector(4)
    n_der = Matrix(2,4)
    shapes.Compute(curve_geometry.Knots, t)

    for i in range(shapes.NbNonzeroPoles):
        n_0[i] = shapes(0, i)
        n_der[0,i] = shapes(1, i)
        n_der[1,i] = shapes(2, i)

    node_indices = [curve.Node(j).Id for j in
        range(shapes.FirstNonzeroPole, shapes.LastNonzeroPole + 1)]


    # Generierung der Elemente pro Integrationspunkt
    element = model_part.CreateNewElement('IgaBeamElement', n+1, node_indices, element_properties)
    element.SetValue(INTEGRATION_WEIGHT, weight)
    element.SetValue(SHAPE_FUNCTION_VALUES, n_0)
    element.SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, n_der)

# Freiheitsgrade einfügen
dof_node = 4
VariableUtils().AddDof(DISPLACEMENT_X, REACTION_X, model_part)
VariableUtils().AddDof(DISPLACEMENT_Y, REACTION_Y, model_part)
VariableUtils().AddDof(DISPLACEMENT_Z, REACTION_Z, model_part)
VariableUtils().AddDof(DISPLACEMENT_ROTATION, REACTION_ROTATION, model_part)

# Randbedingungen: Auflager
# Kontrollpunkt 1
node_1.Fix(DISPLACEMENT_X)
node_1.Fix(DISPLACEMENT_Y)
node_1.Fix(DISPLACEMENT_Z)
node_1.Fix(DISPLACEMENT_ROTATION)

# Randbedingeungen: Knotenlast
load_properties = model_part.GetProperties()[2]
#                             typ,                     Id,  Knoten                   , Eigenschaften
model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [node_4.Id], load_properties)


# Löser konfigurieren
model_part.SetBufferSize(1)

# Verfaren
time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

linear_solver = new_linear_solver_factory.ConstructSolver(Parameters(
    r'{"solver_type": "SkylineLUFactorizationSolver"}'))

# abbruchkriterium
relative_tolerance = 1e-7
absolute_tolerance = 1e-7
conv_criteria = ResidualCriteria(relative_tolerance, absolute_tolerance)
conv_criteria.SetEchoLevel(2)

# Löser
maximum_iterations = 1000 #!! Wenn der Löser nur eine Iteration durchführt erhälst du eine lineare Lösung > Iterationszahl erhöhen!
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
solver.SetEchoLevel(3)


print("\n\n### Berechnung ###")
print("initiale Verformung")
print("Knoten Z: " + str(model_part.GetNode(4).Z))
print("Knoten Z0: " + str((model_part.GetNode(4).Z0)))


for i in range(100):
    F = i * 500/100
    # node_2.SetSolutionStepValue(POINT_LOAD_Y, 1000 * (i + 1) / 10)
    model_part.GetNode(4).SetSolutionStepValue(POINT_LOAD_Z, F)

    # aktuellen modellzustand kopieren
    model_part.CloneTimeStep(i+1)

    # aktuellen zustand lösen
    solver.Solve()

    print("\n\nStep " + str(i) + " :: Knoten Z: " + str(model_part.GetNode(4).Z))
    print("Dehnung: " + str(model_part.GetNode(4).Z - model_part.GetNode(4).Z0))
    print("Verschiebung Z an Knoten 4: " + str(molde_part.GetNode(4).GetValue(DISPLACEMENT_Z)))






print("Prozes time:  %s seconds ---" % (time.time() - start_time))
print('done!')
