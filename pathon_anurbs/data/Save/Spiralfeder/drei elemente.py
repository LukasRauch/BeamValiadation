from kratos_utilities import Model
import os
import ANurbsDev as an
import numpy as np

print('Process ID: ', os.getpid())
print(' ')

geometry = an.Model()
geometry.Load(r'data/spiralfeder.iga')
open('OutputAD.txt', 'w').close()

model = Model(geometry)

for curve_geometry_ptr in geometry.GetByType('CurveGeometry3D'):
    model.add_beam(curve_geometry_ptr)

youngs_modulus = 2e8
Iy = 8.33e-7
Iz = 8.33e-9
model.add_beam_properties('material_1',
    area = 1e-3,   # Querschnittsfläche
    it = 3.0e-8,
    iy = Iz,   # Flächenträgheitsmoment Iy
    iz = Iy,   # Flächenträgheitsmoment Iz
    youngs_modulus = youngs_modulus,
    shear_modulus = youngs_modulus / (2+2*0.00),
)

beam_a = model.beam('curve_a')
beam_b = model.beam('curve_b')
beam_c = model.beam('curve_c')

beam_a.add_stiffness(
    material='material_1',
)
beam_b.add_stiffness(
    material='material_1',
)
beam_c.add_stiffness(
    material='material_1',
)
penalty_dict = {"displacement_x" : 1e+10,
                "displacement_y" : 1e+10,
                "displacement_z" : 1e+10, 
                "torsion" : 1e+10,
                "rotation_2" : 1e+10,
                "rotation_3" : 1e+10}

beam_a.add_coupling(t=beam_a.t1(), other=beam_b, other_t=beam_b.t0(), penalty = penalty_dict, geometry=geometry)
beam_b.add_coupling(t=beam_b.t1(), other=beam_c, other_t=beam_c.t0(), penalty = penalty_dict, geometry=geometry)


beam_a.add_support(t = beam_a.t0(), penalty = {"displacement_x" : 1e+09,
                                                "displacement_y" : 1e+09,
                                                "displacement_z" : 1e+09, 
                                                "torsion" : 1e+09,
                                                "rotation_2" : 1e+09,
                                                "rotation_3" : 1e+09})

# beam_b.add_support(t = beam_b.t0(), penalty = {"displacement_x" : 0,
#                                                 "displacement_y" : 0,
#                                                 "displacement_z" : 0, 
#                                                 "torsion" : 1e09,
#                                                 "rotation_2" : 0,
#                                                 "rotation_3" : 0})

# from KratosMultiphysics import *
# from KratosMultiphysics.IgaApplication import *
# beam_a.nodes[0].Fix(DISPLACEMENT_X)
# beam_a.nodes[0].Fix(DISPLACEMENT_Y)
# beam_a.nodes[0].Fix(DISPLACEMENT_Z)
# beam_a.nodes[0].Fix(DISPLACEMENT_ROTATION)
# beam_a.nodes[1].Fix(DISPLACEMENT_X)
# beam_a.nodes[1].Fix(DISPLACEMENT_Y)
# beam_a.nodes[1].Fix(DISPLACEMENT_Z)
# beam_a.nodes[1].Fix(DISPLACEMENT_ROTATION)


# beam_b.add_node_load(index=-1)

beam_c.add_node_load_moment(t=beam_c.t1(), force=[0, 0, 0], moment=[0,0,1])

beam_a.evaluate_point(t=beam_a.t0(),material='material_1')
beam_a.evaluate_point(t=beam_a.t1(),material='material_1')
beam_b.evaluate_point(t=beam_b.t0(),material='material_1')
beam_b.evaluate_point(t=beam_b.t1(),material='material_1')
beam_c.evaluate_point(t=beam_c.t0(),material='material_1')
beam_c.evaluate_point(t=beam_c.t1(),material='material_1')


model.init_solver()

for step, lam in enumerate(np.linspace(0, 1, 37)[1:]):
    beam_a.make_header(step)

    print(' ##############################')
    print('Force lambda-z: ', lam )  
    M_lam = (2*np.pi*youngs_modulus*Iy/15) * lam
    print('Load = ',  M_lam)

    # beam_b.set_node_load(index=-1, load=[0, 0, 0])


    model.solve(-M_lam)

    beam_a.write_displacement(step)
    beam_b.write_displacement(step)
    beam_c.write_displacement(step)
    # beam_b.write_displacement(step)
    beam_a.write_results(step)
    geometry.Save(r'data/output.iga')


# --- output

beam_a.print_forces(1)
print('BEAM_A')
beam_a.print_displacement()

geometry.Save(r'data/output.iga')

print('Done!')
