from kratos_utilities import Model
import os
import ANurbsDev as an
import numpy as np

print('Process ID: ', os.getpid())
print(' ')

geometry = an.Model()
geometry.Load(r'data/model_2curve_torque.iga')
open('OutputAD.txt', 'w').close()

model = Model(geometry)

for curve_geometry_ptr in geometry.GetByType('CurveGeometry3D'):
    model.add_beam(curve_geometry_ptr)

a = 0.1       # Querschnittshöhe
b = 0.01       # Querschnittsbreite
c1 = 1/3 * ( 1 - 0.63/(a/b) + 0.052/((a/b)**5))
youngs_modulus = 200000
model.add_beam_properties('material_1',
    area = 0.01,   # Querschnittsfläche
    # it = c1*a*b**3,
    it = 0.1,
    iy = 0.0001,   # Flächenträgheitsmoment Iy
    iz = 0.01,   # Flächenträgheitsmoment Iz
    youngs_modulus = youngs_modulus,
    shear_modulus = youngs_modulus / (2+2*0.00),
)

# youngs_modulus = 70000
# model.add_beam_properties('material_2',
#     area = 0.314159,   # Querschnittsfläche
#     it = 5.0e-07,
#     iy = 5e-010,   # Flächenträgheitsmoment Iy
#     iz = 5.0e-06,   # Flächenträgheitsmoment Iz
#     youngs_modulus = youngs_modulus,
#     shear_modulus = youngs_modulus / 2,
# )


beam_a = model.beam('curve_a')
beam_b = model.beam('curve_b')

beam_a.add_stiffness(
    material='material_1',
)

beam_b.add_stiffness(
    material='material_1',
)

penalty_dict = {"displacement_x" : 1e+10,
                "displacement_y" : 1e+10,
                "displacement_z" : 1e+10, 
                "torsion" : 1e+10,
                "rotation_2" : 1e+10,
                "rotation_3" : 1e+10}

beam_a.add_coupling(t=beam_a.t1(), other=beam_b, other_t=beam_b.t0(), penalty = penalty_dict, geometry=geometry)

beam_a.add_support(t = beam_a.t0(), penalty = {"displacement_x" : 1e+10,
                                                "displacement_y" : 1e+10,
                                                "displacement_z" : 1e+10, 
                                                "torsion" : 1e+10,
                                                "rotation_2" : 1e+10,
                                                "rotation_3" : 1e+10})

# beam_b.add_support(t = beam_b.t1(), penalty = {"displacement_x" : 1e+10,
#                                                 "displacement_y" : 1e+10,
#                                                 "displacement_z" : 1e+10, 
#                                                 "torsion" : 1e+10,
#                                                 "rotation_2" : 1e+10,
#                                                 "rotation_3" : 1e+10})

# from KratosMultiphysics import *
# from KratosMultiphysics.IgaApplication import *
# beam_a.nodes[0].Fix(DISPLACEMENT_X)
# beam_a.nodes[0].Fix(DISPLACEMENT_Y)
# beam_a.nodes[0].Fix(DISPLACEMENT_Z)
# beam_a.nodes[0].Fix(DISPLACEMENT_ROTATION)
# beam_a.nodes[1].Fix(DISPLACEMENT_X)
# beam_a.nodes[1].Fix(DISPLACEMENT_Y)
# beam_a.nodes[1].Fix(DISPLACEMENT_Z)
# # beam_a.nodes[1].Fix(DISPLACEMENT_ROTATION)

beam_a.add_node_load(index=-1)
beam_b.add_node_load(index=-1)

# beam_b.add_node_load_moment(t=beam_b.t1(), force=[0, 0, 0], moment=[0,0,0])

# beam_a.evaluate_point(t=beam_a.t0(),material='material_1')
# beam_a.evaluate_point(t=beam_a.t1(),material='material_1')

model.init_solver()

for step, lam in enumerate(np.linspace(0, 1, 51)[1:]):
    beam_a.make_header(step)

    print(' \n##############################')
    print('Force lambda: ', lam )
    F = lam
    M_lam = (2*np.pi/10) * lam
    print('Moment = ',  M_lam)

    # beam_a.set_node_value(index=0, directions=['rotation'], value=F)

    # beam_a.set_node_load(index=-1, load=[ -F ,0, 0])
    beam_a.set_node_load(index=-1, load=[0, 0, -F])
    # beam_b.set_node_load(index=-1, load=[0, F, 0])
    # w = (F * 4**3)/(48 * youngs_modulus * 1)     # Einfeldträger
    # w = (F*10)**2/(youngs_modulus * 1)          # Kragarmpaper
    # print('w = ' + '%.12f' % (w) )      

    # beam_b.set_node_load(index=-1, load=[0,0 , -lam * 0.001])

    model.solve()

    # beam_a.write_displacement(step)
    # beam_b.write_displacement(step)
    # beam_a.write_results(step)
    geometry.Save(r'data/output.iga')


# --- output

beam_a.print_forces(1)
beam_a.print_displacement()

geometry.Save(r'data/output.iga')

print('Done!')
