# module importieren
from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import new_linear_solver_factory

# modell erzeugen
model_part = ModelPart('Model')

# variablen definieren die jeder knoten speichern soll
model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
model_part.AddNodalSolutionStepVariable(REACTION)
model_part.AddNodalSolutionStepVariable(POINT_LOAD)

# elementeigenschaften definieren
truss_properties = model_part.GetProperties()[1] # property-id = 1
truss_properties.SetValue(CROSS_AREA      , 0.01  )
truss_properties.SetValue(PRESTRESS_CAUCHY, 0     )
truss_properties.SetValue(YOUNG_MODULUS   , 210000)
truss_properties.SetValue(POISSON_RATIO   , 0     )
truss_properties.SetValue(DENSITY         , 7856  )

# geometrie erstellen

# knoten                          id  x    y    z
node_1 = model_part.CreateNewNode( 1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode( 2, 0.5, 0.0, -0.5)
node_3 = model_part.CreateNewNode( 3, 1.0, 0.0, 0.0)
# node_4 = model_part.CreateNewNode( 4, 1.0, 0.0, 0.0)


# nurbs gewichte setzen
node_1.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
node_2.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
node_3.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
# node_4.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)

# curve
curve = NodeCurveGeometry3D(
    Degree=1,
    NumberOfNodes=3,
)

curve.SetKnot(Index=0, Value=0.0)
curve.SetKnot(Index=1, Value=1.0)
curve.SetKnot(Index=2, Value=2.0)
# curve.SetKnot(Index=3, Value=3.0)

curve.SetNode(Index=0, Value=node_1)
curve.SetNode(Index=1, Value=node_2)
curve.SetNode(Index=2, Value=node_3)
# curve.SetNode(Index=3, Value=node_4)

# integrationspunkte ermitteln
integration_points = [integration_point for span in curve.Spans()
    for integration_point in IntegrationPoints.Points1D(curve.Degree + 1,
    span)]

# hilfsobjekt um formfunktionen zu berechnen
shapes = CurveShapeEvaluator(Degree=curve.Degree, Order=1)

for i, (t, weight) in enumerate(integration_points):
    # formfunktionen berechnen
    shapes.Compute(curve.Knots(), t)

    # ids der knoten die einen einfluss an der stelle t haben
    node_indices = [curve.Node(j).Id for j in
        range(shapes.FirstNonzeroPole, shapes.LastNonzeroPole + 1)]

    #                                     elementtyp         id     knoten-ids    eigenschaften
    element = model_part.CreateNewElement('IgaTrussElement', i + 1, node_indices, truss_properties)

    # formfunktion aus shape-objekt in vektor kopieren
    n_0 = Vector(shapes.NumberOfNonzeroPoles)
    for i in range(shapes.NumberOfNonzeroPoles):
        n_0[i] = shapes(0, i)

    # ableitung der formfunktion in matrix kopieren
    n_1 = Matrix(shapes.NumberOfNonzeroPoles, 1)
    for i in range(shapes.NumberOfNonzeroPoles):
        n_1[i, 0] = shapes(1, i)

    # integrationsdaten am element speichern
    element.SetValue(INTEGRATION_WEIGHT, weight)
    element.SetValue(SHAPE_FUNCTION_VALUES, n_0)
    element.SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, n_1)

# freiheitsgrade einfügen
VariableUtils().AddDof(DISPLACEMENT_X, REACTION_X, model_part)
VariableUtils().AddDof(DISPLACEMENT_Y, REACTION_Y, model_part)
VariableUtils().AddDof(DISPLACEMENT_Z, REACTION_Z, model_part)

# randbedingungen: auflager
node_1.Fix(DISPLACEMENT_X)
node_1.Fix(DISPLACEMENT_Y)
node_1.Fix(DISPLACEMENT_Z)

node_2.Fix(DISPLACEMENT_Y)
# node_2.Fix(DISPLACEMENT_Z)

node_3.Fix(DISPLACEMENT_X)
node_3.Fix(DISPLACEMENT_Y)
node_3.Fix(DISPLACEMENT_Z)

# randbedingungen: knotenlasten
load_properties = model_part.GetProperties()[2]

#                             typ                       id knoten       eigenschaften
model_part.CreateNewCondition('PointLoadCondition3D1N', 1, [node_2.Id], load_properties)


# löser konfigurieren

# anzahl der schritte die gespeichert werden sollen
model_part.SetBufferSize(1)

# verfahren
time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

linear_solver = new_linear_solver_factory.ConstructSolver(Parameters(
    r'{"solver_type": "SkylineLUFactorizationSolver"}'))


# abbruchkriterium
relative_tolerance = 1e-7
absolute_tolerance = 1e-7
conv_criteria = ResidualCriteria(relative_tolerance, absolute_tolerance)
conv_criteria.SetEchoLevel(0)

# löser
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


# system lösen

f = open("Outputfile.md","w")
f.write("\n\n### Berechnung ###")
f.write("initiale Verformung")
f.write("\nKnoten Z: " + str(node_2.Z))
f.write("\nKnoten Z0: " + str(node_2.Z0))

# Print to Ecxel
import xlsxwriter

workbook = xlsxwriter.Workbook('OutputData.xlsx')
worksheet = workbook.add_worksheet()

bold = workbook.add_format({'bold': True})
num_fmt = workbook.add_format({'num_format': '$0,00000'})
# num_fmt.set_num_format({'#,##0'})


worksheet.write('A1', 'Python Output', bold)
worksheet.write('A2', 'KRATOS engin', bold)
worksheet.write('A4', 'Knot Nr.')
# worksheet.write('C4', 'Pointload.')
# worksheet.write('B4', 'Position Z')



for i in range(100):
    F = i * 500/100
    # node_2.SetSolutionStepValue(POINT_LOAD_Y, 1000 * (i + 1) / 10)
    node_2.SetSolutionStepValue(POINT_LOAD_Z, F)

    # aktuellen modellzustand kopieren
    model_part.CloneTimeStep(i+1)

    # aktuellen zustand lösen
    solver.Solve()
    
    f.write("\n\nStep " + str(i) + " :: Knoten Z: " + str(node_2.Z))
    f.write("\nDehnung: " + str(node_2.Z - node_2.Z0))
    # f.write("\nVerschiebung X an Knoten 2: " + str(node_2.GetValue(Z)))

    worksheet.write('A' + str(5+i) , str(i)  ) 
    worksheet.write('C' + str(5+i) , float(F), num_fmt ) 
    worksheet.write('B' + str(5+i) , float(node_2.Z) , num_fmt)


workbook.close()

# auflagerverschiebung am knoten 2
print(node_2.Z - node_2.Z0)


# open("Outputfile.md")

