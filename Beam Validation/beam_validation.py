# module importieren
from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import new_linear_solver_factory
import os
import yaml
import json
import time

start_time = time.time()

print('Process ID:', os.getpid())

stream = open(r'C:/_Masterarbeit/Debug/Outputdata.yaml', 'r')


for data in yaml.load_all(stream):
    print('\n% ')
    print('Prozess Step #', data['t'])
    print('//...')
    # modell erzeugen
    model_part = ModelPart('Model')

    # variablen definieren die jeder knoten speichern soll
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_ROTATION)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(REACTION_ROTATION)
    model_part.AddNodalSolutionStepVariable(POINT_LOAD)

    # elementeigenschaften definieren
    element_properties = model_part.GetProperties()[1] # property-id = 1
    element_properties.SetValue(CROSS_AREA, data['section_area'])
    element_properties.SetValue(YOUNG_MODULUS, data['material_youngsmodulus'])
    element_properties.SetValue(SHEAR_MODULUS, data['material_shearmodulus'])
    element_properties.SetValue(MOMENT_OF_INERTIA_Y, data["section_momentofinertia_y"])
    element_properties.SetValue(MOMENT_OF_INERTIA_Z, data["section_momentofinertia_z"])
    element_properties.SetValue(MOMENT_OF_INERTIA_T, data["section_momentofinertia_t"])
      
    coords = data['coordinates']
    displacements = data['displacements']
    nNode_Size= data['nNode_Size']
    node_indices = []
    nodes = []
    
    kratos_curve = NodeCurveGeometry3D(Degree = 3, NumberOfNodes = 4)

    # Zuweisung der Knotenkoordinaten     
    for i in range(nNode_Size):         # Knotenkoordintaten initial Displacement
                                        #  Id X                Y                Z 
        nodes = model_part.CreateNewNode( i+1, coords[i][0] , coords[i][1] , coords[i][2] )
        node_indices.append(nodes.Id)  
        
    # for i in range(nNode_Size):
    # for i, (u, v, w, r) in enumerate(data['displacements']):
        nodes.SetSolutionStepValue(DISPLACEMENT_X, displacements[i][0])
        nodes.SetSolutionStepValue(DISPLACEMENT_Y, displacements[i][1])
        nodes.SetSolutionStepValue(DISPLACEMENT_Z, displacements[i][2])
        nodes.SetSolutionStepValue(DISPLACEMENT_ROTATION, displacements[i][3])


###################
## TODO ControlPoint Weight einfügen
## TODO Degree festelegen 
###################

    # # Erzeugung der n-1 Elemente [Elemententyp , Id ,[ Node_Ids ],  Element Eigenschaften]
    # for i in range( nNode_Size - 1 ):
    #     node_indices = list(range(1, nNode_Size-1))
    element = model_part.CreateNewElement('IgaBeamElement', 1, node_indices, element_properties)

    # Externe Elementvorgaben
    element.SetValue(T0, data['T0_vec'])
    element.SetValue(N0, data['N0_vec'])
    element.SetValue(PHI, data['PHI_PROP'])
    element.SetValue(PHI_0_DER, data['PHI_DER_PROP'])
    element.SetValue(INTEGRATION_WEIGHT, data['gp_w']*10/2) # *2  
    # element.SetValue(KNOTSPANN, data['kNodspann'])

    # Shape Functions
    n_0 = Vector(4)
    n_deriv = Matrix(2,4)
    # n_2 = Matrix(4,1)
    shape_0 = data['shape_0']
    shape_1 = data['shape_1']
    shape_2 = data['shape_2']
    
    for i in range(4): 
        n_0[i] = shape_0[i]
        n_deriv[0,i] = shape_1[i]
        n_deriv[1,i] = shape_2[i]
        # n_2[i,0] = shape_2[i]
        
    element.SetValue(SHAPE_FUNCTION_VALUES, n_0)
    element.SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, n_deriv) 
    # element.SetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES, n_2)
  
    
    element.SetValue(DEBUG_EXPECTED_DATA, json.dumps(data))

    # freiheitsgrade einfügen
    
## TODO Belegung der voreingestllen freiheitsgrade 
    dof_node = 4

    VariableUtils().AddDof(DISPLACEMENT_X, REACTION_X, model_part)
    VariableUtils().AddDof(DISPLACEMENT_Y, REACTION_Y, model_part)
    VariableUtils().AddDof(DISPLACEMENT_Z, REACTION_Z, model_part)
    if(dof_node >= 4): 
        VariableUtils().AddDof(DISPLACEMENT_ROTATION, REACTION_ROTATION, model_part)

    

#########################
## TODO Randbedingungen einfügen
#########################
    nodes.Fix(DISPLACEMENT_X)
    nodes.Fix(DISPLACEMENT_Y)
    nodes.Fix(DISPLACEMENT_Z)
    if(dof_node >= 4): nodes.Fix(DISPLACEMENT_ROTATION)

    
############################
##! TODO Lasten einfügen
############################
    # Dummy Last
    i = 4
    nodes.SetSolutionStepValue(POINT_LOAD_X, 1000 * (1) / 10)

    
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
    maximum_iterations = 0
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

    model_part.CloneTimeStep(0)

    solver.Solve()

print('finished!')
print("Prozes time:  %s seconds ---" % (time.time() - start_time))