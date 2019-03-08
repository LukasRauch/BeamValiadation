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
from scipy import integrate
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

# model = an.Model.open(r'C:\_Masterarbeit\BeamValidation\Balken\Balken.iga')
model = an.Model.open(r'C:\_Masterarbeit\BeamValidation\Balken\Balken_schief.iga')

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
phi = 0                         # manuelle Vorgabe der Rotation
phi_der = 0                     # manuelle Vorgabe der Rotation 1st Ableitung
knots = curve_geometry.Knots

element_count = 0

def calculate_tau(var_t):
    """returns tau depending on the parameter t"""
    _, r_1, r_2, r_3 = kratos_curve.DerivativesAt(T = var_t , Order = 3)
    a = np.cross(r_1,r_2)
    d = np.dot(np.cross(r_1, r_2), r_3)
    delta = np.linalg.norm(r_1)
    alpha = np.linalg.norm(a)
    kappa = alpha / (delta **(3))
    tau   = d / (alpha ** 2)
    return tau

def calculate_delta(var_t):
    """returns delta depending on the parameter t"""
    _, r_1, r_2, r_3 = kratos_curve.DerivativesAt(T = var_t , Order = 3)
    a = np.cross(r_1,r_2)
    d = np.dot(np.cross(r_1, r_2), r_3)
    delta = np.linalg.norm(r_1)
    return delta

function = lambda var_t: calculate_tau(var_t) * calculate_delta(var_t)

def normalize(v):
    return v / np.linalg.norm(v)

def normalize_1(v, v_1):
    return (-(np.dot(v, v_1) + np.dot(v_1, v)) * v + 2 * np.dot(v, v) * v_1) / (2 * np.dot(v, v)**(3/2))

def cross_1(u, u_1, v, v_1):
    return np.cross(u, v_1) + np.cross(u_1, v)

for n, (t, weight) in enumerate(integration_points):    # 11 Integrationspunkte
    element_count += 1
    # Geometie
    theta_0 = 0*np.pi/180
    _, r_1, r_2, r_3 = kratos_curve.DerivativesAt(T = t , Order = 3)

    def base_vector():
        tol = 1e-10
        a = np.cross(r_1,r_2)
        a_1 = cross_1(r_1, r_2, r_2, r_3)                
        d = np.dot(a, r_3)     
        delta = np.linalg.norm(r_1)             
        alpha = np.linalg.norm(a)  
        if delta >= tol:             
            kappa = alpha / delta**3
            tau   = d / alpha**2
            # tau   = np.linalg.det([r_1,r_2,r_3]) / alpha**2

            T = r_1 / delta     # Tangente A1 
            B = a / alpha       # Transversale A3
            N = np.cross(B, T)  # Normale A2

            T_1 = normalize_1(r_1, r_2)
            B_1 = normalize_1(a, a_1)
            N_1 = cross_1(B, B_1, T, T_1)

            # theta integration

            def func(t):
                _, r_1, r_2, r_3 = kratos_curve.DerivativesAt(T = t , Order = 3)
                a = np.cross(r_1, r_2)   
                d = np.dot(a, r_3)
                delta = np.linalg.norm(r_1)             
                alpha = np.linalg.norm(a)               
                tau = d / alpha**2
                return -tau * delta

            # theta = theta_0 +  tau * delta * weight 
            # theta_1 =  tau * delta 
            theta   = theta_0 + integrate.romberg(func, 0, t)
            theta_1 = tau * delta

            A1 = r_1
            A1_1 = r_2

            A2 = (  N * np.cos(theta) + B * np.sin(theta)).tolist()         # check sollte so weit funktionieren    
            A2_1 = ( ( + N_1 * np.cos(theta) - N * np.sin(theta) * theta_1
                    + B_1 * np.sin(theta) + B * np.cos(theta) * theta_1)).tolist()
            
            A3 = ( -N * np.sin(theta) + B * np.cos(theta)).tolist()         # check sollte so weit funktionieren
            A3_1 =   ( - N_1 * np.sin(theta) - N * np.cos(theta) * theta_1
                        + B_1 * np.cos(theta) - B * np.sin(theta) * theta_1).tolist()
            passtol = 1e-10
        a = np.cross(r_1,r_2)
        a_1 = cross_1(r_1, r_2, r_2, r_3)                
        d = np.dot(a, r_3)     
        delta = np.linalg.norm(r_1)             
        alpha = np.linalg.norm(a)  
        if alpha >= tol:             
            kappa = alpha / delta**3
            tau   = d / alpha**2
            # tau   = np.linalg.det([r_1,r_2,r_3]) / alpha**2

            T = r_1 / delta     # Tangente A1 
            B = a / alpha       # Transversale A3
            N = np.cross(B, T)  # Normale A2

            T_1 = normalize_1(r_1, r_2)
            B_1 = normalize_1(a, a_1)
            N_1 = cross_1(B, B_1, T, T_1)

            # theta integration

            def func(t):
                _, r_1, r_2, r_3 = kratos_curve.DerivativesAt(T = t , Order = 3)
                a = np.cross(r_1, r_2)   
                d = np.dot(a, r_3)
                delta = np.linalg.norm(r_1)             
                alpha = np.linalg.norm(a)               
                tau = d / alpha**2
                return -tau * delta

            # theta = theta_0 +  tau * delta * weight 
            # theta_1 =  tau * delta 
            theta   = theta_0 + integrate.romberg(func, 0, t)
            theta_1 = tau * delta

            A1 = r_1
            A1_1 = r_2

            A2 = (  N * np.cos(theta) + B * np.sin(theta)).tolist()         # check sollte so weit funktionieren    
            A2_1 = ( ( + N_1 * np.cos(theta) - N * np.sin(theta) * theta_1
                    + B_1 * np.sin(theta) + B * np.cos(theta) * theta_1)).tolist()
            
            A3 = ( -N * np.sin(theta) + B * np.cos(theta)).tolist()         # check sollte so weit funktionieren
            A3_1 =   ( - N_1 * np.sin(theta) - N * np.cos(theta) * theta_1
                        + B_1 * np.cos(theta) - B * np.sin(theta) * theta_1).tolist()
            pass
        else:   # Fall: Gerader Stab:: Krümmung = inf
            T = r_1 / delta     # Tangente A1
            A1 = r_1
            A1_1 = r_2
            if np.array_equal(T,[1,0,0]):
                A2 = [0,1,0]
            elif np.array_equal(T,[0,1,0]):
                A2 = [0,0,1]
            elif np.array_equal(T,[0,0,1]):
                A2 = [1,0,0]
            else:
                A2   = (normalize([-r_1[1] , r_1[0] , 0])).tolist()

            A3   = (normalize(np.cross(A2,T))).tolist()

        A2_1 = [0,0,0]
        A3_1 = [0,0,0]

        return A1, A1_1, A2, A2_1, A3, A3_1


    A1, A1_1, A2, A2_1, A3, A3_1 = base_vector()


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

    # Generierung der Elemente pro Integrationspunkt
    # element = model_part.CreateNewElement('IgaBeamElement', n+1, node_indices, element_properties)
    element = model_part.CreateNewElement('IgaBeamADElement', n+1, node_indices, element_properties)
    element.SetValue(INTEGRATION_WEIGHT, weight)
    element.SetValue(SHAPE_FUNCTION_VALUES              , n_0)     # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DER_1         , n_1)     # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DER_2         , n_2)     # Typ Vektor
    element.SetValue(SHAPE_FUNCTION_LOCAL_DER_3         , n_3)     # Typ Vektor
    element.SetValue(T0                                 , r_1)
    element.SetValue(T0_DER                             , r_2)

    element.SetValue(BASE_A1                            , A1)
    element.SetValue(BASE_A2                            , A2)
    element.SetValue(BASE_A3                            , A3)
    element.SetValue(BASE_A1_1                          , A1_1)
    element.SetValue(BASE_A2_1                          , A2_1)
    element.SetValue(BASE_A3_1                          , A3_1)
    ### manuelle Vor
    element.SetValue(N0                                 , n0)
    element.SetValue(PHI                                , phi)
    element.SetValue(PHI_DER_1                          , phi_der)


# Randbedingungen: Knotenlast
load_properties = model_part.GetProperties()[2] # propperty-ID = 2
#                             typ,                     Id,  Knoten                   , Eigenschaften
model_part.CreateNewCondition('PointLoadCondition3D1N', 2, [model_part.GetNode(curve_geometry.NbPoles).Id], load_properties)
# _________________________________________________________________________________________________________________
# # # # # Definition: Bettung
# position_t = 0


# # Formfunktionen an Gausspunkt n
# n_0 = Vector(4)                                     
# n_1 = Vector(4)                                     
# n_2 = Vector(4)                                     
# n_3 = Vector(4)                                     
# # n_der = Matrix(2,4)                                 
# shapes.Compute(curve_geometry.Knots, position_t)
# node_indices = np.arange(shapes.FirstNonzeroPole+1, shapes.LastNonzeroPole+2, dtype=int)

# theta_0 = 0*np.pi/180
# point, r_1, r_2, r_3 = kratos_curve.DerivativesAt(T = t , Order = 3)

# A1, A1_1, A2, A2_1, A3, A3_1 = base_vector()

# for i in range(shapes.NbNonzeroPoles):
#     n_0[i] = shapes(0, i)
#     n_1[i] = shapes(1, i)
#     n_2[i] = shapes(2, i)
#     n_3[i] = shapes(3, i)

# # # Tangentenvektor ausgewertet an Gausspunkt n
# # # Normierung des Tangentenvektors erfolgt Kratos-intern
# point, r_1, r_2 = kratos_curve.DerivativesAt(T=position_t, Order=2)  # Tangentenvektor am aktuellen Integrationspunkt auswerten

# # # Generierung der Elemente pro Integrationspunkt
# # # element = model_part.CreateNewElement('IgaBeamADElement', n+1, node_indices, element_properties)
# element_dirichlet_condition = model_part.CreateNewElement('IgaBeamWeakDirichletCondition', element_count+1, node_indices, element_properties)
# element_dirichlet_condition.SetValue(INTEGRATION_WEIGHT                 , 1)  # *2
# element_dirichlet_condition.SetValue(SHAPE_FUNCTION_VALUES              , n_0)     # Typ Vektor
# element_dirichlet_condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_1         , n_1)     # Typ Vektor
# element_dirichlet_condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_2         , n_2)     # Typ Vektor
# element_dirichlet_condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_3         , n_3)     # Typ Vektor
# element_dirichlet_condition.SetValue(T0                                 , r_1)
# element_dirichlet_condition.SetValue(T0_DER                             , r_2)

# element_dirichlet_condition.SetValue(BASE_A1                            , A1)
# element_dirichlet_condition.SetValue(BASE_A2                            , A2)
# element_dirichlet_condition.SetValue(BASE_A3                            , A3)
# element_dirichlet_condition.SetValue(BASE_A1_1                          , A1_1)
# element_dirichlet_condition.SetValue(BASE_A2_1                          , A2_1)
# element_dirichlet_condition.SetValue(BASE_A3_1                          , A3_1)

# ### manuelle Vorgabe
# element_dirichlet_condition.SetValue(N0                                 , n0)
# element_dirichlet_condition.SetValue(PHI                                , phi)
# element_dirichlet_condition.SetValue(PHI_DER_1                          , phi_der)
# ### Randbedingungen 
# element_dirichlet_condition.SetValue(PENALTY_DISPLACEMENT               , 1e12)
# element_dirichlet_condition.SetValue(PENALTY_ROTATION                   , 1e12)
# element_dirichlet_condition.SetValue(PENALTY_TORSION                    , 1e12)
# element_dirichlet_condition.SetValue(DIRICHLET_CONDITION_TYPE           , 123)    # 1 Displacement, 2 Torsion , 3 Rotation Winkel, 4 Steigung

# # # # # _________________________________________________________________________________________________________________
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

# # model_part.GetNode(2).Fix(DISPLACEMENT_X)
model_part.GetNode(2).Fix(DISPLACEMENT_Y)
model_part.GetNode(2).Fix(DISPLACEMENT_Z)
# model_part.GetNode(1).Fix(DISPLACEMENT_ROTATION)



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

num_pole = curve_geometry.NbPoles
num_load_steps = 10

disp_X = []
disp_Y = []
disp_Z = []

disp_X = np.empty([num_load_steps+1, num_pole])
disp_Y = np.empty([num_load_steps+1, num_pole])
disp_Z = np.empty([num_load_steps+1, num_pole])

# PLOT _________________________________________________________________________________
multi_curve = Multi.MultiCurve()


for i in range(0, num_load_steps+1):
    F           = i * 15/num_load_steps
    # moment_vec  = [i * 10/num_load_steps ,0 , 0]

    model_part.GetNode(curve_geometry.NbPoles).SetSolutionStepValue(POINT_LOAD_Z, F)
    # model_part.GetElement(n+2)
    # element_load_properties.SetValue(LOAD_VECTOR_MOMENT, moment_vec)

    # aktuellen modellzustand kopieren
    model_part.CloneTimeStep(i+1)

    # aktuellen zustand lösen
    print("solver step: ", i)
    solver.Solve()


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

    # # Draw the control points polygon, the 3D curve and the vectors
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

# Set evaluation delta
multi_curve.delta = 0.001
# multi_curve.evaluate
# plot the controlpoint polygon and the evaluated curve
vis_comp = VisMPL.VisCurve3D()
multi_curve.vis = vis_comp
multi_curve
multi_curve.render(cpcolor='black', evalcolor='red')
# multi_curve.render()