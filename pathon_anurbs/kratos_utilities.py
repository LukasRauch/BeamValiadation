from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import new_linear_solver_factory
import ANurbsDev as an
import numpy as np
from scipy import integrate
import yaml
import json
import pandas as pd


def normalized(v):
    return v / np.linalg.norm(v)

def normalized_1(v, v_1):
    return v_1 / np.linalg.norm(v) - np.dot(v, v_1) * v / np.linalg.norm(v)**3

def cross_1(u, u_1, v, v_1):
    return np.cross(u, v_1) + np.cross(u_1, v)

def cross_v_identity(v):
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])

def compute_rod(v, phi):
    return np.cos(phi) * np.eye(3) + cross_v_identity(np.sin(phi) * v)

def compute_rod_1(v, v_1, phi, phi_1):
    return np.cos(phi) * phi_1 * cross_v_identity(v) + np.sin(phi) * (cross_v_identity(v_1) - phi_1 * np.eye(3))

def compute_lam(v1, v2):
    v1_x_v2 = np.cross(v1, v2)
    v1_d_v2 = np.dot(v1, v2)

    lam = v1_d_v2 * np.eye(3) + cross_v_identity(v1_x_v2)

    if v1_d_v2 + 1.0 > 1e-7:
        lam = lam + np.outer(v1_x_v2, v1_x_v2) * 1.0 / (1.0 + v1_d_v2)
    else:
        l_v1_x_v2 = np.linalg.norm(v1_x_v2)

        e_hat = v1_x_v2

        if l_v1_x_v2 > 1e-7:
            e_hat = e_hat / l_v1_x_v2
            lam = lam + np.outer(e_hat, e_hat) * (1 - v1_d_v2)

    return lam

def compute_lam_1(v1, v1_1, v2, v2_1):
    T0_T   = np.dot(v1, v2)
    T0_T_1 = np.dot(v1, v2_1) + np.dot(v1_1, v2)
    T0xT   = np.cross(v1, v2)
    T0xT_1 = np.cross(v1, v2_1) + np.cross(v1_1, v2)

    d = 1.0 / (1.0 + T0_T)

    o = np.outer(T0xT_1, T0xT) + np.outer(T0xT, T0xT_1)

    return T0_T_1 * np.eye(3) + cross_v_identity(T0xT_1) - T0_T_1 * np.power(d, 2) * np.outer(T0xT, T0xT) + d * o


class Model:
    def __init__(self, geometry):
        # modell erzeugen
        model_part = ModelPart('Model')

        # variablen definieren die jeder knoten speichern soll
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT_ROTATION)
        model_part.AddNodalSolutionStepVariable(REACTION)
        model_part.AddNodalSolutionStepVariable(REACTION_ROTATION)
        model_part.AddNodalSolutionStepVariable(POINT_LOAD)

        self.model_part = model_part

        self.nodes = []
        self._node_keys = {}

        self.beams = []
        self._beam_keys = {}

        self.properties = []
        self._property_keys = {}

        self.elements = []
        self._element_keys = {}

        self._time_step = 0

        self.geometry = geometry

        self.penalty = []
        self.penalty_key = {}

    def add_node(self, key, location):
        node_id = len(self.nodes) + 1
        node = self.model_part.CreateNewNode(node_id, *location)
        node.AddDof(DISPLACEMENT_X, REACTION_X)
        node.AddDof(DISPLACEMENT_Y, REACTION_Y)
        node.AddDof(DISPLACEMENT_Z, REACTION_Z)
        node.AddDof(DISPLACEMENT_ROTATION, REACTION_ROTATION)
        self.nodes.append(node)
        self._node_keys[key] = node
        return node

    def node(self, key):
        return self._node_keys[key]

    def add_beam(self, curve_geometry_ptr):
        beam = Beam(self, curve_geometry_ptr)
        self.beams.append(beam)
        self._beam_keys[beam.key] = beam

    def beam(self, key):
        return self._beam_keys[f'{key}.CurveGeometry3D']

    def add_properties(self, key=None):
        property_id = len(self.properties) + 1

        property = self.model_part.GetProperties()[property_id]

        self.properties.append(property)

        if key is not None:
            self._property_keys[key] = property

    def add_beam_properties(self, key, area, it, iy, iz, youngs_modulus, shear_modulus):
        property_id = len(self.properties) + 1

        property = self.model_part.GetProperties()[property_id]
        property.SetValue(CROSS_AREA, area)               # m²
        property.SetValue(MOMENT_OF_INERTIA_T, it)        # m4
        property.SetValue(MOMENT_OF_INERTIA_Y, iy)        # m4
        property.SetValue(MOMENT_OF_INERTIA_Z, iz)        # m4
        property.SetValue(YOUNG_MODULUS, youngs_modulus)  # kN/m²
        property.SetValue(SHEAR_MODULUS, shear_modulus)   # kN/m²

        self.properties.append(property)
        self._property_keys[key] = property

    def property(self, key):
        return self._property_keys[key]

    def add_element(self, element_type, nodes, property):
        element_id = len(self.elements) + 1
        node_ids = [node.Id for node in nodes]
        element = self.model_part.CreateNewElement(element_type, element_id, node_ids, property)
        self.elements.append(element)
        self._element_keys[element_id] = element
        return element

    def add_condition(self, condition_type, nodes, property):
        element_id = len(self.elements) + 1
        node_ids = [node.Id for node in nodes]
        element = self.model_part.CreateNewCondition(condition_type, element_id, node_ids, property)
        self.elements.append(element)
        self._element_keys[element_id] = element
        return element

    def solve(self):
        # Löser konfigurieren
        self.model_part.SetBufferSize(1)

        # Verfahren
        time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        linear_solver = new_linear_solver_factory.ConstructSolver(Parameters(
            r'{"solver_type": "eigen_sparse_lu"}'))

        # Abbruchkriterium
        relative_tolerance = 1e-7
        absolute_tolerance = 1e-7
        # conv_criteria = ResidualCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria = DisplacementCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria.SetEchoLevel(1)

        # Löser
        maximum_iterations = 300 #!! Wenn der Löser nur eine Iteration durchführt erhälst du eine lineare Lösung > Iterationszahl erhöhen!
        compute_reactions = True
        reform_dofs_at_each_iteration = True
        move_mesh_flag = True

        solver = ResidualBasedNewtonRaphsonStrategy(
            self.model_part,
            time_scheme,
            linear_solver,
            conv_criteria,
            maximum_iterations,
            compute_reactions,
            reform_dofs_at_each_iteration,
            move_mesh_flag
        )
        solver.SetEchoLevel(1)

        self._time_step += 1

        self.model_part.CloneTimeStep(self._time_step)

        solver.Solve()

        self.update()

    def update(self):
        for beam in self.beams:
            beam.update(self._time_step)

class Beam:
    def __init__(self, model, curve_geometry_ptr):

        self.clear_memory()

        curve_geometry = curve_geometry_ptr.Data()

        nodes = []

        for i, pole in enumerate(curve_geometry.Poles()):
            node_key = (curve_geometry_ptr.Key(), i)

            node = model.add_node(node_key, location=pole)

            nodes.append(node)

        self.model = model
        self.curve_geometry_ptr = curve_geometry_ptr
        self.curve_geometry = curve_geometry
        self.nodes = nodes

    @property
    def key(self):
        return self.curve_geometry_ptr.Key()

    @property
    def model_part(self):
        return self.model.model_part

    @property
    def t0(self):
        return self.curve_geometry.Domain().T0()

    def fix_node(self, index, directions):
        node = self.nodes[index]

        if 'x' in directions:
            node.Fix(DISPLACEMENT_X)
        if 'y' in directions:
            node.Fix(DISPLACEMENT_Y)
        if 'z' in directions:
            node.Fix(DISPLACEMENT_Z)
        if 'rotation' in directions:
            node.Fix(DISPLACEMENT_ROTATION)

    def add_node_load(self, index, load=[0,0,0]):
        node = self.nodes[index]

        load_properties = self.model.add_properties()

        self.model.add_condition('PointLoadCondition3D1N', [node], load_properties)

        node.SetSolutionStepValue(POINT_LOAD_X, load[0])
        node.SetSolutionStepValue(POINT_LOAD_Y, load[1])
        node.SetSolutionStepValue(POINT_LOAD_Z, load[2])

    def set_node_load(self, index, load):
        node = self.nodes[index]

        node.SetSolutionStepValue(POINT_LOAD_X, load[0])
        node.SetSolutionStepValue(POINT_LOAD_Y, load[1])
        node.SetSolutionStepValue(POINT_LOAD_Z, load[2])

        #Print load vector

        geometry = self.model.geometry

        scale = 1
        if np.amax(np.absolute(load)) != 0:
            scale = np.amax(np.absolute(load))

        load_0 = [node.X, node.Y, node.Z]
        load_vec = [node.X - load[0]//scale, node.Y - load[1]/scale, node.Z - load[2]/scale]

        line_ptr = geometry.Add(an.Line3D(a=load_vec, b=load_0))
        line_ptr.Attributes().SetLayer(f'LoadVector')
        line_ptr.Attributes().SetColor(f'#ff0000')
        line_ptr.Attributes().SetArrowhead('End')

    def _func(self, t):
        curve_geometry = self.curve_geometry
        _, r_1, r_2, r_3 = curve_geometry.DerivativesAt(t = t , order = 3)
        a = np.cross(r_1, r_2)
        d = np.dot(a, r_3)
        delta = np.linalg.norm(r_1)
        alpha = np.linalg.norm(a)
        tau = d / alpha**2
        return -tau * delta

    def frame_at(self, t):
        curve_geometry = self.curve_geometry
        p, r_1, r_2, r_3 = curve_geometry.DerivativesAt(t = t, order = 3)

        # print("p: " , p)

        tol = 1e-10
        a = np.cross(r_1,r_2)
        a_1 = cross_1(r_1, r_2, r_2, r_3)
        d = np.dot(a, r_3)
        delta = np.linalg.norm(r_1)
        alpha = np.linalg.norm(a)
        if alpha >= tol:
            kappa = alpha / delta**3
            tau   = d / alpha**2

            T = r_1 / delta     # Tangente a1
            B = a / alpha       # Transversale a3
            N = np.cross(B, T)  # Normale a2

            T_1 = normalized_1(r_1, r_2)
            B_1 = normalized_1(a, a_1)
            N_1 = cross_1(B, B_1, T, T_1)

            theta   = integrate.romberg(self._func, 0, t,divmax=10)
            theta_1 = tau * delta

            a1 = r_1
            a1_1 = r_2

            a2 = (  N * np.cos(theta) + B * np.sin(theta))        
            a2_1 = ( ( + N_1 * np.cos(theta) - N * np.sin(theta) * theta_1
                    + B_1 * np.sin(theta) + B * np.cos(theta) * theta_1))

            a3 = ( -N * np.sin(theta) + B * np.cos(theta))        
            a3_1 =   ( - N_1 * np.sin(theta) - N * np.cos(theta) * theta_1
                        + B_1 * np.cos(theta) - B * np.sin(theta) * theta_1)
            pass
        else:   # Fall: Gerader Stab:: Krümmung = inf
            T = r_1 / delta    
            a1 = r_1
            a1_1 = r_2
            if np.array_equal(T,[0,0,1]):
                a2 = np.array([0,1,0])
            else:
                a2   = normalized([-r_1[1] , r_1[0] , 0])

            a3   = np.cross(T,a2)

        a2_1 = np.array([0,0,0])
        a3_1 = np.array([0,0,0])

        return p, a1, a1_1, a2, a2_1, a3, a3_1

    def add_moment(self, t, vector, material):
        if not isinstance(property, Properties):
            material = self.model.property(material)

        curve_geometry = self.curve_geometry
        # model_part = self.model_part

        # integration_degree = curve_geometry.Degree() + 1

        nonzero_node_indices, shape_functions = curve_geometry.ShapeFunctionsAt(t, order=3)

        nonzero_nodes = [self.nodes[index] for index in nonzero_node_indices]

        condition = self.model.add_element('IgaBeamMomentCondition', nonzero_nodes, material)

        condition.SetValue(INTEGRATION_WEIGHT, 1)     # FIXME: integration_weight von der Kratoskurve beziehen.

        condition.SetValue(SHAPE_FUNCTION_VALUES     , shape_functions[0])
        condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_1, shape_functions[1])
        condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_2, shape_functions[2])
        condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_3, shape_functions[3])

        _, A1, A1_1, A2, A2_1, A3, A3_1 = self.frame_at(t)

        condition.SetValue(BASE_A1, A1.tolist())
        condition.SetValue(BASE_A2, A2.tolist())
        condition.SetValue(BASE_A3, A3.tolist())
        condition.SetValue(BASE_A1_1, A1_1.tolist())
        condition.SetValue(BASE_A2_1, A2_1.tolist())
        condition.SetValue(BASE_A3_1, A3_1.tolist())

        condition.SetValue(LOAD_VECTOR_MOMENT, vector)


    def add_support(self, t, penalty):

        geometry = self.model.geometry

        material = self.model.add_beam_properties('dummy_material',
            area = 0, it = 0, iy = 0, iz = 0,
            youngs_modulus = 0, shear_modulus = 0,
        )

        bool_support_x = False
        bool_support_y = False
        bool_support_z = False

        curve_geometry = self.curve_geometry
        model_part = self.model_part

        integration_degree = curve_geometry.Degree() + 1

        curve_geometry = self.curve_geometry

        nonzero_node_indices, shape_functions = curve_geometry.ShapeFunctionsAt(t, order=3)

        nonzero_nodes = [self.nodes[index] for index in nonzero_node_indices]

        condition = self.model.add_element('IgaBeamWeakDirichletCondition', nonzero_nodes, material)
        
        condition.SetValue(SHAPE_FUNCTION_VALUES     , shape_functions[0])
        condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_1, shape_functions[1])
        condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_2, shape_functions[2])
        condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_3, shape_functions[3])

        p, A1, A1_1, A2, A2_1, A3, A3_1 = self.frame_at(t)

        condition.SetValue(BASE_A1, A1.tolist())
        condition.SetValue(BASE_A2, A2.tolist())
        condition.SetValue(BASE_A3, A3.tolist())
        condition.SetValue(BASE_A1_1, A1_1.tolist())
        condition.SetValue(BASE_A2_1, A2_1.tolist())
        condition.SetValue(BASE_A3_1, A3_1.tolist())

        DISPLACEMENT_X = 0    # default
        DISPLACEMENT_Y = 0    # default
        DISPLACEMENT_Z = 0    # default
        TORSION = 0         # default
        ROTATION = 0        # default

        if 'displacement_x' in penalty:
            DISPLACEMENT_X = penalty["displacement_x"]
            bool_support_x = True
        if 'disp_x' in penalty:
            DISPLACEMENT_X = penalty["disp_x"]
            bool_support_x = True
        if 'displacement_y' in penalty:
            DISPLACEMENT_Y = penalty["displacement_y"]
            bool_support_y = True
        if 'disp_y' in penalty:
            DISPLACEMENT_Y = penalty["disp_y"]
            bool_support_y = True
        if 'displacement_z' in penalty:
            DISPLACEMENT_Z = penalty["displacement_z"]
            bool_support_z = True
        if 'disp_z' in penalty:
            DISPLACEMENT_Z = penalty["disp_z"]
            bool_support_z = True
        if 'torsion' in penalty:
            TORSION = penalty["torsion"]
        if 'tors' in penalty:
            TORSION = penalty["tors"]
        if 'rotation_2' in penalty:
            ROTATION_2 = penalty["rotation_2"]
        if 'rot_2' in penalty:
            ROTATION_2 = penalty["rot_2"]
        if 'rotation_3' in penalty:
            ROTATION_3 = penalty["rotation_3"]
        if 'rot_3' in penalty:
            ROTATION_3 = penalty["rot_3"]

        condition.SetValue(PENALTY_DISPLACEMENT_X, DISPLACEMENT_X)
        condition.SetValue(PENALTY_DISPLACEMENT_Y, DISPLACEMENT_Y)
        condition.SetValue(PENALTY_DISPLACEMENT_Z, DISPLACEMENT_Z)
        condition.SetValue(PENALTY_TORSION, TORSION)
        condition.SetValue(PENALTY_ROTATION_2, ROTATION_2)
        condition.SetValue(PENALTY_ROTATION_3, ROTATION_3)

        if bool_support_x:
            line_ptr = geometry.Add(an.Line3D(a=np.add(p,np.array([0.01,0,0])), b=p))
            line_ptr.Attributes().SetLayer(f'Support')
            line_ptr.Attributes().SetColor(f'#00ff00')
            line_ptr.Attributes().SetArrowhead('End')

        
        if bool_support_y:
            line_ptr = geometry.Add(an.Line3D(a=np.add(p,np.array([0,0.01,0])), b=p))
            line_ptr.Attributes().SetLayer(f'Support')
            line_ptr.Attributes().SetColor(f'#00ff00')
            line_ptr.Attributes().SetArrowhead('End')

    
        if bool_support_z:
            line_ptr = geometry.Add(an.Line3D(a=np.add(p,np.array([0,0,0.01])), b=p))
            line_ptr.Attributes().SetLayer(f'Support')
            line_ptr.Attributes().SetColor(f'#00ff00')
            line_ptr.Attributes().SetArrowhead('End')


    def add_stiffness(self, material):
        if not isinstance(property, Properties):
            material = self.model.property(material)

        curve_geometry = self.curve_geometry
        # model_part = self.model_part

        integration_points = []

        integration_degree = curve_geometry.Degree() + 1

        for span in curve_geometry.Spans():
            if span.Length() < 1e-7:
                continue

            integration_points += an.IntegrationPoints.Points1D(
                degree=integration_degree,
                domain=span,
            )

        for i, (t, weight) in enumerate(integration_points):
            nonzero_node_indices, shape_functions = curve_geometry.ShapeFunctionsAt(t, order=3)

            nonzero_nodes = [self.nodes[index] for index in nonzero_node_indices]

            element = self.model.add_element('IgaBeamADElement', nonzero_nodes, material)

            element.SetValue(INTEGRATION_WEIGHT, weight)

            element.SetValue(SHAPE_FUNCTION_VALUES     , shape_functions[0])
            element.SetValue(SHAPE_FUNCTION_LOCAL_DER_1, shape_functions[1])
            element.SetValue(SHAPE_FUNCTION_LOCAL_DER_2, shape_functions[2])
            element.SetValue(SHAPE_FUNCTION_LOCAL_DER_3, shape_functions[3])

            p, A1, A1_1, A2, A2_1, A3, A3_1 = self.frame_at(t)

            element.SetValue(BASE_A1, A1.tolist())
            element.SetValue(BASE_A2, A2.tolist())
            element.SetValue(BASE_A3, A3.tolist())
            element.SetValue(BASE_A1_1, A1_1.tolist())
            element.SetValue(BASE_A2_1, A2_1.tolist())
            element.SetValue(BASE_A3_1, A3_1.tolist())

            element.SetValue(GAUSS_POINT, list(p))

            frame = open('frames.txt', 'a')
            frame.write( str(p.tolist()) + str(A1.tolist()) + str(A2.tolist()) + str(A3.tolist()) + '\n')
            frame.close()

    def t0(self):
        act_curve_geometry  = self.curve_geometry.Clone()
        domain = act_curve_geometry.Domain()
        return domain.T0()

    def t1(self):
        act_curve_geometry  = self.curve_geometry.Clone()
        domain = act_curve_geometry.Domain()
        return domain.T1()
 

    def add_coupling(self, t, other, other_t, penalty, geometry):
        material = self.model.add_beam_properties('dummy_material',
            area = 0, it = 0, iy = 0, iz = 0,
            youngs_modulus = 0, shear_modulus = 0, 
        )

        curve_geometry_a = self.curve_geometry
        curve_geometry_b = other.curve_geometry
        
        nonzero_node_indices_a, shape_functions_a = curve_geometry_a.ShapeFunctionsAt(t = t, order=3)
        nonzero_node_indices_b, shape_functions_b = curve_geometry_b.ShapeFunctionsAt(other_t, order=3)

        nonzero_nodes = []

        nonzero_nodes = [self.nodes[index] for index in nonzero_node_indices_a]
        nonzero_nodes_b = [other.nodes[index] for index in nonzero_node_indices_b]

        nonzero_nodes.extend(nonzero_nodes_b)

        condition = self.model.add_element('IgaBeamADWeakCoupling', nonzero_nodes, material)

        condition.SetValue(SHAPE_FUNCTION_VALUES     , shape_functions_a[0])
        condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_1, shape_functions_a[1])
        condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_2, shape_functions_a[2])
        condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_3, shape_functions_a[3])

        condition.SetValue(SHAPE_FUNCTION_VALUES_B     , shape_functions_b[0])
        condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_1_B, shape_functions_b[1])
        condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_2_B, shape_functions_b[2])
        condition.SetValue(SHAPE_FUNCTION_LOCAL_DER_3_B, shape_functions_b[3])

        XA, A1, A1_1, A2, A2_1, A3, A3_1 = self.frame_at(t)

        condition.SetValue(BASE_A1, A1.tolist())
        condition.SetValue(BASE_A2, A2.tolist())
        condition.SetValue(BASE_A3, A3.tolist())
        condition.SetValue(BASE_A1_1, A1_1.tolist())
        condition.SetValue(BASE_A2_1, A2_1.tolist())
        condition.SetValue(BASE_A3_1, A3_1.tolist())

        XB, B1, B1_1, B2, B2_1, B3, B3_1 = other.frame_at(other_t)

        condition.SetValue(BASE_B1, B1.tolist())
        condition.SetValue(BASE_B2, B2.tolist())
        condition.SetValue(BASE_B3, B3.tolist())
        condition.SetValue(BASE_B1_1, B1_1.tolist())

        geometry.Add(an.Point3D(location=XA))
        geometry.Add(an.Point3D(location=XB))

        if 'displacement_x' in penalty:
            DISPLACEMENT_X = penalty["displacement_x"]
        if 'disp_x' in penalty:
            DISPLACEMENT_X = penalty["disp_x"]
        if 'displacement_y' in penalty:
            DISPLACEMENT_Y = penalty["displacement_y"]
        if 'disp_y' in penalty:
            DISPLACEMENT_Y = penalty["disp_y"]
        if 'displacement_z' in penalty:
            DISPLACEMENT_Z = penalty["displacement_z"]
        if 'disp_z' in penalty:
            DISPLACEMENT_Z = penalty["disp_z"]
        if 'torsion' in penalty:
            TORSION = penalty["torsion"]
        if 'tors' in penalty:
            TORSION = penalty["tors"]
        if 'rotation_2' in penalty:
            ROTATION_2 = penalty["rotation_2"]
        if 'rot_2' in penalty:
            ROTATION_2 = penalty["rot_2"]
        if 'rotation_3' in penalty:
            ROTATION_3 = penalty["rotation_3"]
        if 'rot_3' in penalty:
            ROTATION_3 = penalty["rot_3"]

        condition.SetValue(PENALTY_DISPLACEMENT_X, DISPLACEMENT_X)
        condition.SetValue(PENALTY_DISPLACEMENT_Y, DISPLACEMENT_Y)
        condition.SetValue(PENALTY_DISPLACEMENT_Z, DISPLACEMENT_Z)
        condition.SetValue(PENALTY_TORSION, TORSION)
        condition.SetValue(PENALTY_ROTATION_2, ROTATION_2)
        condition.SetValue(PENALTY_ROTATION_3, ROTATION_3)

    def update(self, time_step):
        geometry = self.model.geometry
        act_curve_geometry = self.curve_geometry.Clone()

        for i, node in enumerate(self.nodes):
            act_curve_geometry.SetPole(i, [node.X, node.Y, node.Z])

        act_curve_geometry_key = f'{self.key}.Step<{time_step}>'
        act_curve_geometry_ptr = geometry.Add(act_curve_geometry_key, act_curve_geometry)

        key = self.key.replace('.CurveGeometry3D', '')
        curve_key = f'{key}.Step<{time_step}>'
        curve_ptr = geometry.Add(curve_key, an.Curve3D(act_curve_geometry_ptr))

        curve_ptr.Attributes().SetLayer(f'Step<{time_step}>')

        domain = act_curve_geometry.Domain()

        for t in np.linspace(domain.T0(), domain.T1(), 10):
            _, A1, A1_1, A2, A2_1, A3, A3_1 = self.frame_at(t)

            T = normalized(A1)
            T_1 = normalized_1(A1, A1_1)

            def act_phi_at(t):
                nonzero_node_indices, shape_functions = self.curve_geometry.ShapeFunctionsAt(t, order=1)

                nonzero_nodes = [self.nodes[index] for index in nonzero_node_indices]

                act_values = np.array([node.GetSolutionStepValue(DISPLACEMENT_ROTATION) for node in nonzero_nodes])

                phi, phi_1 = np.dot(shape_functions, act_values)

                return phi, phi_1

            phi, phi_1 = act_phi_at(t)

            x, a1, a1_1 = act_curve_geometry.DerivativesAt(t, order=2)

            a11 = np.dot(a1, a1)
            a = np.sqrt(a11)

            t = normalized(a1)
            t_1 = a1_1 / a - np.dot(a1, a1_1) * a1 / a**3

            rod = compute_rod(t, phi)
            rod_1 = compute_rod_1(t, t_1, phi, phi_1)

            lam = compute_lam(T, t)
            lam_1 = compute_lam_1(T, T_1, t, t_1)

            rod_lam = np.dot(rod  , lam  )
            rod_1_lam = np.dot(rod_1, lam  )
            rod_lam_1 = np.dot(rod  , lam_1)

            a2 = np.dot(rod_lam, A2)
            a3 = np.dot(rod_lam, A3)

            scale = 0.25

            # print('a2', a2, A2)
            # print('a3', a3, A3)

            # line_ptr = geometry.Add(an.Line3D(a=x, b=x+a1*scale/np.linalg.norm(a1)))
            # line_ptr.Attributes().SetLayer(f'Step<{time_step}>')
            # line_ptr.Attributes().SetColor(f'#ff0000')
            # # line_ptr.Attributes().SetArrowhead('End')

            line_ptr = geometry.Add(an.Line3D(a=x, b=x+a2*scale))
            line_ptr.Attributes().SetLayer(f'Step<{time_step}>')
            line_ptr.Attributes().SetColor(f'#00ff00')
            # line_ptr.Attributes().SetArrowhead('End')

            line_ptr = geometry.Add(an.Line3D(a=x, b=x+a3*scale))
            line_ptr.Attributes().SetLayer(f'Step<{time_step}>')
            line_ptr.Attributes().SetColor(f'#0000ff')
            # line_ptr.Attributes().SetArrowhead('End')

    def print_forces(self, scale):
        fname = 'cutting_force.txt'
        data = np.loadtxt(fname, dtype={'names': ('Id', 'x', 'y', 'z', 'N', 'M2', 'M3', 't_0', 't_1', 't_2', 'a2_0', 'a2_1', 'a2_2', 'a3_0', 'a3_1', 'a3_2'), 
                                        'formats': ('i4', 'f4', 'f4' , 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')})

        frame = np.loadtxt('frames.txt', dtype=np.str, delimiter='\s')
        
        geometry = self.model.geometry

        i = len(data) -1
        list_n  = np.array([])
        list_m2 = np.array([])
        list_m3 = np.array([])

        while data[i][0] >= 0:
            list_n = np.append(list_n , np.absolute(data[i][4]))
            list_m2 = np.append(list_m2, np.absolute(data[i][5]))
            list_m3 = np.append(list_m3, np.absolute(data[i][6]))

            if data[i][0] == 1: 
                break
            i -= 1
        
        norm_n =np.amax(list_n)
        norm_m2 = np.amax(list_m2)
        norm_m3 = np.amax(list_m3)

        a2 = [data['a2_0'][-1], data['a2_1'][-1], data['a2_2'][-1]]        
        a3 = [data['a3_0'][-1], data['a3_1'][-1], data['a3_2'][-1]]

        n  = np.dot(a3, data['N'][-1] * scale)
        if norm_n != 0:  n  = np.dot(a3, data['N'][-1] * scale/norm_n) 

        m2 = np.dot(a2, data['M2'][-1] * scale)
        if norm_m2 != 0: m2 = np.dot(a2, data['M2'][-1] * scale / norm_m2)

        m3 = np.dot(a3, data['M3'][-1] * scale)
        if norm_m3 != 0: m3 = np.dot(a3, data['M3'][-1] * scale / norm_m3)

        #Print starting line off the loop
        x_old_n  = np.add([data[-1][1], data[-1][2], data[-1][3]], n)
        x_old_m2 = np.add([data[-1][1], data[-1][2], data[-1][3]], m2)
        x_old_m3 = np.add([data[-1][1], data[-1][2], data[-1][3]], m3) 

        line_ptr = geometry.Add(an.Line3D(a=np.add(x_old_n,-n), b=np.add(np.add(x_old_n,-n), n)))
        line_ptr.Attributes().SetLayer(f'Normalkraft N')
        if data[i][4] <= 0:
            line_ptr.Attributes().SetColor(f'#ff0000')     # negative forces = red   
        else:
            line_ptr.Attributes().SetColor(f'#0000ff')     # negative forces = red    

        line_ptr = geometry.Add(an.Line3D(a=np.add(x_old_m2, -m2), b=np.add(np.add(x_old_m2, -m2), m2)))
        line_ptr.Attributes().SetLayer(f'Moment Mz')
        if data[i][5] <= 0:
            line_ptr.Attributes().SetColor(f'#ff0000')
        else:
            line_ptr.Attributes().SetColor(f'#0000ff')

        line_ptr = geometry.Add(an.Line3D(a=np.add(x_old_m3, -m3), b=np.add(np.add(x_old_m3, -m3), m3)))
        line_ptr.Attributes().SetLayer(f'Moment My')
        if data[i][6] <= 0:
            line_ptr.Attributes().SetColor(f'#ff0000')
        else:
            line_ptr.Attributes().SetColor(f'#0000ff')


        i = len(data) -2
        while data[i][0] >= 0:

            a2 = [data['a2_0'][i], data['a2_1'][i], data['a2_2'][i]]
            a3 = [data['a3_0'][i], data['a3_1'][i], data['a3_2'][i]]

            x = [data[i][1], data[i][2], data[i][3]]

            n  = np.dot(a3, data['N'][i] * scale)
            if norm_n != 0:  n  = np.dot(a3, data['N'][i] * scale/norm_n) 

            m2 = np.dot(a2, data['M2'][i] * scale)
            if norm_m2 != 0: m2 = np.dot(a2, data['M2'][i] * scale / norm_m2)

            m3 = np.dot(a3, data['M3'][i] * scale)
            if norm_m3 != 0: m3 = np.dot(a3, data['M3'][i] * scale / norm_m3)

            line_ptr = geometry.Add(an.Line3D(a=x, b=np.add(x, n)))
            line_ptr.Attributes().SetLayer(f'Normalkraft N')
            if data[i][4] <= 0:
                line_ptr.Attributes().SetColor(f'#ff0000')
            else:
                line_ptr.Attributes().SetColor(f'#0000ff')

            line_ptr = geometry.Add(an.Line3D(a=x_old_n, b=np.add(x, n)))
            line_ptr.Attributes().SetLayer(f'Normalkraft N')
            if data[i][4] <= 0:
                line_ptr.Attributes().SetColor(f'#ff0000')
            else:
                line_ptr.Attributes().SetColor(f'#0000ff')
            x_old_n = np.add(x, n)

            line_ptr = geometry.Add(an.Line3D(a=x, b=np.add(x, m2)))
            line_ptr.Attributes().SetLayer(f'Moment Mz')
            if data[i][5] <= 0:
                line_ptr.Attributes().SetColor(f'#ff0000')
            else:
                line_ptr.Attributes().SetColor(f'#0000ff')

            line_ptr = geometry.Add(an.Line3D(a=x_old_m2, b=np.add(x, m2)))
            line_ptr.Attributes().SetLayer(f'Moment Mz')
            if data[i][5] <= 0:
                line_ptr.Attributes().SetColor(f'#ff0000')
            else:
                line_ptr.Attributes().SetColor(f'#0000ff')
            x_old_m2 = np.add(x, m2)

            line_ptr = geometry.Add(an.Line3D(a=x, b=np.add(x, m3)))
            line_ptr.Attributes().SetLayer(f'Moment My')
            if data[i][6] <= 0:
                line_ptr.Attributes().SetColor(f'#ff0000')
            else:
                line_ptr.Attributes().SetColor(f'#0000ff')

            line_ptr = geometry.Add(an.Line3D(a=x_old_m3, b=np.add(x, m3)))
            line_ptr.Attributes().SetLayer(f'Moment My')
            if data[i][6] <= 0:
                line_ptr.Attributes().SetColor(f'#ff0000')
            else:
                line_ptr.Attributes().SetColor(f'#0000ff')
            x_old_m3 = np.add(x, m3)
            

            if data[i][0] == 1:
                break

            i -= 1

    def print_displacement(self, Id=[]):
        act_curve_geometry = self.curve_geometry.Clone()

        #Header
        print('\nDisplacements of ' + str(self.key))

        if Id:
            node = self.nodes[Id]
            print(
                f"{'Id:':<4}{k:<4}",
                f"{'x:':<4}{node.X - node.X0:<30}" ,
                f"{'y:':<4}{node.Y - node.Y0:<30}" ,
                f"{'y:':<4}{node.Z - node.Z0:<30}" 
                )

        else:
            for k, pole in enumerate(act_curve_geometry.Poles()):
                node = self.nodes[k]
                print(
                    f"{'Id:':<4}{k:<4}",
                    f"{'x:':<4}{node.X - node.X0:<30}" ,
                    f"{'y:':<4}{node.Y - node.Y0:<30}" ,
                    f"{'y:':<4}{node.Z - node.Z0:<30}" 
                    )

    def write_displacement(self):
        act_curve_geometry = self.curve_geometry.Clone()

        #Header
        data = open('displacements.txt', 'a')
        data.write('Displacements of ' + str(self.key) + '\n')
        data.close()

        for k, pole in enumerate(act_curve_geometry.Poles()):
            node = self.nodes[k]
            displacement = {k: {'Id': k+1,
                                'x': node.X - node.X0,
                                'y': node.Y - node.Y0,
                                'z': node.Z - node.Z0, 
                                }
                            }
                                    # node = self.nodes[index]

            df = pd.DataFrame.from_dict([displacement])
            df.to_csv('displacements.txt', header=False, index=False, mode='a')

        data = open('displacements.txt', 'a')
        data.write('\n\n')
        data.close()

    def clear_memory(self):
        open('cutting_force.txt', 'w').close()
        open('frames.txt', 'w').close()
        open('displacements.txt', 'w')


