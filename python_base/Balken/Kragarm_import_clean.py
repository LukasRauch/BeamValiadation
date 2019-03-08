# IMPORT MODUL
from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import new_linear_solver_factory
import ANurbs as an

# PROCESS INFO 
print('Process ID: ', os.getpid())
print(' ')


#__________________________________________________________________________________________________________
# IMPORT RHINO GEOMETRY
model = an.Model.open(r'C:\_Masterarbeit\BeamValidation\Balken\Balken.iga')

curve_item = model.of_type('Curve3D')[0]
curve = curve_item.data
curve_geometry = curve_item.geometry().geometry

# GENERATE MODEL
model_part = ModelPart('Model')

# DEFINE ELEMENT VARIABLES
model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
model_part.AddNodalSolutionStepVariable(DISPLACEMENT_ROTATION)
model_part.AddNodalSolutionStepVariable(REACTION)
model_part.AddNodalSolutionStepVariable(REACTION_ROTATION)
model_part.AddNodalSolutionStepVariable(POINT_LOAD) 

# DEFINE ELEMENT PROPPERTIES
element_properties = model_part.GetProperties()[1]     # elementpropertie - ID '1'
element_properties.SetValue(CROSS_AREA         ,   100)
element_properties.SetValue(YOUNG_MODULUS      ,   1)
element_properties.SetValue(SHEAR_MODULUS      ,   0.5)
element_properties.SetValue(MOMENT_OF_INERTIA_Y,   100)
element_properties.SetValue(MOMENT_OF_INERTIA_Z,   500)
element_properties.SetValue(MOMENT_OF_INERTIA_T,   100)
element_properties.SetValue(POISSON_RATIO      ,   0)
element_properties.SetValue(DENSITY             ,   78.5)
element_properties.SetValue(POISSON_RATIO       ,   0)

# DEFINE CURVE 
curve = NodeCurveGeometry3D(Degree = curve_geometry.Degree , NumberOfNodes = curve_geometry.NbPoles )

# GENERATE GEOMETRY
node_indices = []
for i in range(curve_geometry.NbPoles):
    # Kontrollpunkte                ID , X                         , Y                         , Z
    node = model_part.CreateNewNode(i+1, curve_geometry.Pole(i)[0],curve_geometry.Pole(i)[1], curve_geometry.Pole(i)[2])
    node_indices.append(i+1)
    # NURBS Gewicht (Hier erst mal noch alle Gewichte 1)
    node.SetValue(NURBS_CONTROL_POINT_WEIGHT, 1.0)
    # Nodes erstellen
    curve.SetNode(Index = i, Value = node)
    
for i, (Knot) in enumerate(curve_geometry.Knots):
    curve.SetKnot(Index = i, Value = Knot)

# REFERENCE GEOMETRY 
n0      = [0, 0, -1]
phi     = 0
phi_der = 0

# COMPUTE SHAPE FUNCTIONS 
integration_points = curve_item.IntegrationPoints()
shapes = an.CurveShapeEvaluator(Degree = curve_geometry.Degree, Order = 2)

for n, (t, weight) in enumerate(integration_points):

    n_0     = Vector(shapes.NbNonzeroPoles, 0)
    n_der   = Matrix(2, shapes.NbNonzeroPoles, 0)
    shapes.Compute(curve_geometry.Knots, t)

    for i in range(shapes.NbNonzeroPoles):
        n_0[i]      = shapes(0, i)
        n_der[0,i]  = shapes(1, i)
        n_der[1,i]  = shapes(2, i)

# GENERATE ELEMENTS 
    point, tangent = curve.DerivativesAt(T=t, Order=1)

    element = model_part.CreateNewElement('IgaBeamElement', n+1, node_indices, element_properties)
    element.SetValue(INTEGRATION_WEIGHT                 , weight)
    element.SetValue(SHAPE_FUNCTION_VALUES              , n_0)
    element.SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES   , n_der)
    element.SetValue(T0                                 , tangent)
    element.SetValue(N0                                 , n0)
    element.SetValue(PHI                                , phi)
    element.SetValue(PHI_0_DER                          , phi_der)

# DEFINE DEGREES OF FREEDOM
VariableUtils().AddDof(DISPLACEMENT_X, REACTION_X, model_part)
VariableUtils().AddDof(DISPLACEMENT_Y, REACTION_Y, model_part)
VariableUtils().AddDof(DISPLACEMENT_Z, REACTION_Z, model_part)
VariableUtils().AddDof(DISPLACEMENT_ROTATION, REACTION_ROTATION, model_part)

# DEFINE BOUNDARY CONDITIONS 
# Einspannung Node 1
model_part.GetNode(1).Fix(DISPLACEMENT_X)
model_part.GetNode(1).Fix(DISPLACEMENT_Y)
model_part.GetNode(1).Fix(DISPLACEMENT_Z)
model_part.GetNode(1).Fix(DISPLACEMENT_ROTATION)
# Verschiebung Node 2
model_part.GetNode(2).Fix(DISPLACEMENT_Y)
model_part.GetNode(2).Fix(DISPLACEMENT_Z)

# DEFINE LOAD PROPERTIES
load_properties = model_part.GetProperties()[2]     # loadpropertie - ID '2'
#                              Typ                     ,ID ,Kontrollpunkt-ID         , LoadProperties
model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [model_part.GetNode(4).Id], load_properties)

# CONFIGURATE SOLVER
# Verfahren
time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

linear_solver = new_linear_solver_factory.ConstructSolver(Parameters(
    # r'{"solver_type": "eigen_sparse_lu"}'))
    r'{"solver_type": "SkylineLUFactorizationSolver"}'))

# Buffer Size
model_part.SetBufferSize(1)

# Abbruchkriterien 
relative_tolerance = 1e-7
absolute_tolerance = 1e-7
conv_criteria = ResidualCriteria(relative_tolerance, absolute_tolerance)

# Echo Level
conv_criteria.SetEchoLevel(1)

# Solver Input
maximum_iterations              = 100
compute_reactions               = True
reform_dofs_at_each_iteration   = True
move_mesh_flag                  = True

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

# SOLVE SYSTEM
# Anzahl Zeitschritte
load_steps = 10
for i in range(load_steps+1):
    # Last Aufbringen
    F = i * 1 / load_steps
    model_part.GetNode(4).SetSolutionStepValue(POINT_LOAD_Z, F)

    # Zeitschritt Klonen
    model_part.CloneTimeStep(i+1)

    # Aktuellen Zeitschritt Lösen
    print("Solver step: ", i)
    solver.Solve()

    # Verschiebung Z Drucken
    print("DSIP-Z = ", model_part.GetNode(4).Z - model_part.GetNode(4).Z0)

print("Exit!")


