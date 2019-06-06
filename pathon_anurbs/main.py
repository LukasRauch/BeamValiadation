from kratos_utilities import Model
import os
import ANurbsDev as an
import numpy as np

# os.environ["OMP_NUM_THREADS"] = "1"
# os.environ["set OMP_NUM_THREADS=1"]


print('Process ID: ', os.getpid())
print(' ')

geometry = an.Model()
geometry.Load(r'data/torsion_arch.iga')
open('OutputAD.txt', 'w').close()

model = Model(geometry)

for curve_geometry_ptr in geometry.GetByType('CurveGeometry3D'):
    model.add_beam(curve_geometry_ptr)

youngs_modulus =1.999e8
Iy = 8.33e-7
Iz = 8.33e-7
It = 8.33e-7
model.add_beam_properties('material_1',
    area = 0.001,   # Querschnittsfläche
    it = It,
    iy = Iy,   # Flächenträgheitsmoment Iy
    iz = Iz,   # Flächenträgheitsmoment Iz
    youngs_modulus = youngs_modulus,
    shear_modulus = youngs_modulus / (2+2*0.25),
)

# model.add_beam_properties('material_2',
#     area = 0.1,   # Querschnittsfläche
#     it = 1e6,
#     iy = 6.0e-1,   # Flächenträgheitsmoment Iy
#     iz = 6.0e-1,   # Flächenträgheitsmoment Iz
#     youngs_modulus = 1e7,
#     shear_modulus = 1e7 / (2+2*0.00),
# )

beam_a = model.beam('curve_a')
# beam_b = model.beam('curve_b')
# beam_c = model.beam('curve_c')

beam_a.add_stiffness(
    material='material_1',
)
# beam_b.add_stiffness(
#     material='material_1',
# )
# beam_c.add_stiffness(
#     material='material_1',
# )

penalty_dict = {"displacement_x" : 1e10,
                "displacement_y" : 1e10,
                "displacement_z" : 1e10, 
                "torsion" : 1e11,
                "rotation_2" : 1e10,
                "rotation_3" : 1e10}

# beam_a.add_coupling(t=beam_a.t1(), other=beam_b, other_t=beam_b.t0(), penalty = penalty_dict, geometry=geometry)
# beam_b.add_coupling(t=beam_b.t1(), other=beam_c, other_t=beam_c.t0(), penalty = penalty_dict, geometry=geometry)

# support = beam_a.add_support(t = beam_a.t0(), penalty = {"displacement_x" : 1e10,
#                                                 "displacement_y" : 1e10,
#                                                 "displacement_z" : 1e10, 
#                                                 "torsion" : 0,
#                                                 "rotation_2" : 0,
#                                                 "rotation_3" : 0})

# support = beam_a.add_support(t = beam_a.t1(), penalty = {"displacement_x" : 1e10,
#                                                 "displacement_y" : 0,
#                                                 "displacement_z" : 1e10, 
#                                                 "torsion" : 0,
#                                                 "rotation_2" : 0,
#                                                 "rotation_3" : 0})


from KratosMultiphysics import *
from KratosMultiphysics.IgaApplication import *
beam_a.nodes[0].Fix(DISPLACEMENT_X)
beam_a.nodes[0].Fix(DISPLACEMENT_Y)
beam_a.nodes[0].Fix(DISPLACEMENT_Z)
beam_a.nodes[0].Fix(DISPLACEMENT_ROTATION)
beam_a.nodes[1].Fix(DISPLACEMENT_X)
beam_a.nodes[1].Fix(DISPLACEMENT_Y)
beam_a.nodes[1].Fix(DISPLACEMENT_Z)
# beam_a.nodes[1].Fix(DISPLACEMENT_ROTATION)


beam_a.add_node_load(index=-1)
# beam_b.add_node_load(index=-1)
# beam_c.add_node_load(index=-1)

# beam_c.add_node_load_moment(t=beam_c.t1(), force=[0, 0, 0], moment=[1,0,0])

# beam_a.evaluate_point(t=beam_a.t0(),material='material_1')
# beam_a.evaluate_point(t=beam_a.t1(),material='material_1')
# beam_b.evaluate_point(t=beam_b.t0(),material='material_2')
# beam_b.evaluate_point(t=beam_b.t1(),material='material_2')
# beam_c.evaluate_point(t=beam_c.t0(),material='material_1')
# beam_c.evaluate_point(t=beam_c.t1(),material='material_1')

model.init_solver()
from KratosMultiphysics.IgaApplication import *


for step, lam in enumerate(np.linspace(0, 1, 41)[1:]):
    beam_a.make_header(step)

    print(' \n##############################')
    print('Force lambda: ', lam )
    F = lam  * np.pi
    M_lam = lam * np.pi * 0.5 * youngs_modulus* It / 8
    print('Moment = ',  M_lam)

    beam_a.set_node_value(index=0, directions=['rotation'], value=F)
    # beam_a.set_node_value(index=1, directions=['rotation'], value=F)

    # A2 = condition.GetValue(BASE_A2)
    # A3 = condition.GetValue(BASE_A3)
    # A2_1 = condition.GetValue(BASE_A2_1) 
    # A3_1 =condition.GetValue(BASE_A3_1) 

    # def rotate(angle, vector):
    #     x = np.cos(angle) * vector[0] - np.sin(angle) * vector[1]
    #     y = np.sin(angle) * vector[0] + np.cos(angle) * vector[1]
    #     z = vector[2]

    #     return [x,y,z]

    # support.SetValue(ALPHA_LOAD, F)


    # beam_a.set_node_load(index=-1, load=[ 0, F, 0])
    # beam_b.set_node_load(index=-1, load=[ 0, 0, F])

    # beam_b.set_node_load(index=-1, load=[0,0 , lam * 0.5 ])

    model.solve()

    beam_a.write_displacement(step)
    # beam_b.write_displacement(step)
    # beam_c.write_displacement(step)
    # beam_b.write_displacement(step)
    # beam_a.write_results(step)l
    geometry.Save(r'data/output.iga')


# --- output

beam_a.print_forces(scale=1)
beam_a.print_displacement()

geometry.Save(r'data/output.iga')

print('Done!')
