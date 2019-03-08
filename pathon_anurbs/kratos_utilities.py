from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import new_linear_solver_factory
import ANurbsDev as an
import numpy as np


def normalized(v):
    return v / np.linalg.norm(v)

def normalized_1(v, v_1):
    return v_1 / np.linalg.norm(v) - np.dot(v, v_1) * v / np.linalg.norm(v)**3

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
        relative_tolerance = 1e-4
        absolute_tolerance = 1e-4
        conv_criteria = DisplacementCriteria(relative_tolerance, absolute_tolerance)
        conv_criteria.SetEchoLevel(1)

        # Löser
        maximum_iterations = 200 #!! Wenn der Löser nur eine Iteration durchführt erhälst du eine lineare Lösung > Iterationszahl erhöhen!
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

    def frame_at(self, t):
        curve_geometry = self.curve_geometry

        # FIXME: für gekrümmte curven -> minimal rotation frame

        p, a1, a1_1 = curve_geometry.DerivativesAt(t, order=2)

        a2 = np.array([0, 1, 0])
        a2_1 = np.array([0, 0, 0])

        a3 = np.array([0, 0, 1])
        a3_1 = np.array([0, 0, 0])

        return p, a1, a1_1, a2, a2_1, a3, a3_1

    def add_stiffness(self, material):
        if not isinstance(property, Properties):
            material = self.model.property(material)

        curve_geometry = self.curve_geometry
        model_part = self.model_part

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

            # Generierung der Elemente pro Integrationspunkt
            # element = model_part.CreateNewElement('IgaBeamElement', n+1, node_indices, element_properties)
            element = self.model.add_element('IgaBeamADElement', nonzero_nodes, material)

            element.SetValue(INTEGRATION_WEIGHT, weight)

            element.SetValue(SHAPE_FUNCTION_VALUES     , shape_functions[0])
            element.SetValue(SHAPE_FUNCTION_LOCAL_DER_1, shape_functions[1])
            element.SetValue(SHAPE_FUNCTION_LOCAL_DER_2, shape_functions[2])
            element.SetValue(SHAPE_FUNCTION_LOCAL_DER_3, shape_functions[3])

            _, A1, A1_1, A2, A2_1, A3, A3_1 = self.frame_at(t)

            element.SetValue(BASE_A1, A1.tolist())
            element.SetValue(BASE_A2, A2.tolist())
            element.SetValue(BASE_A3, A3.tolist())
            element.SetValue(BASE_A1_1, A1_1.tolist())
            element.SetValue(BASE_A2_1, A2_1.tolist())
            element.SetValue(BASE_A3_1, A3_1.tolist())

    def add_support(self, t, penalty):
        pass

    def add_coupling(self, t, other, other_t, penalty):
        curve_geometry_a = self.curve_geometry
        curve_geometry_b = other.curve_geometry
        pass

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

            scale = 0.3

            line_ptr = geometry.Add(an.Line3D(a=x, b=x+a1*scale))
            line_ptr.Attributes().SetLayer(f'Step<{time_step}>')
            line_ptr.Attributes().SetColor(f'#ff0000')
            # line_ptr.Attributes().SetArrowhead('End')

            line_ptr = geometry.Add(an.Line3D(a=x, b=x+a2*scale))
            line_ptr.Attributes().SetLayer(f'Step<{time_step}>')
            line_ptr.Attributes().SetColor(f'#00ff00')
            # line_ptr.Attributes().SetArrowhead('End')
            
            line_ptr = geometry.Add(an.Line3D(a=x, b=x+a3*scale))
            line_ptr.Attributes().SetLayer(f'Step<{time_step}>')
            line_ptr.Attributes().SetColor(f'#0000ff')
            # line_ptr.Attributes().SetArrowhead('End')
