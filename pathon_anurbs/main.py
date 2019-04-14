from kratos_utilities import Model
import os
import ANurbsDev as an
import numpy as np

print('Process ID: ', os.getpid())
print(' ')

geometry = an.Model()
geometry.Load(r'data/einfeldtraeger.iga')
open('OutputAD.txt', 'w').close()

model = Model(geometry)

for curve_geometry_ptr in geometry.GetByType('CurveGeometry3D'):
    model.add_beam(curve_geometry_ptr)

a = 0.01       # Querschnittshöhe
b = 0.01       # Querschnittsbreite
c1 = 1/3 * ( 1 - 0.63/(a/b) + 0.052/((a/b)**5))
youngs_modulus = 100

model.add_beam_properties('material_1',
    area = 200,   # Querschnittsfläche
    it = 0.01,
    iy = 0.1,   # Flächenträgheitsmoment Iy
    iz = 0.1,   # Flächenträgheitsmoment Iz
    youngs_modulus = youngs_modulus,
    shear_modulus = youngs_modulus / 2,
)
model.add_beam_properties('material_2',
    area = 100000,   # Querschnittsfläche
    it = 100,
    iy = 10000,   # Flächenträgheitsmoment Iy
    iz = 10000,   # Flächenträgheitsmoment Iz
    youngs_modulus = youngs_modulus,
    shear_modulus = youngs_modulus / 2,
)


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
                                                "torsion" : 0e+10,
                                                "rotation_2" : 0e+10,
                                                "rotation_3" : 0e+10})

beam_b.add_support(t = beam_b.t1(), penalty = {"displacement_x" : 0e+10,
                                                "displacement_y" : 1e+10,
                                                "displacement_z" : 1e+10, 
                                                "torsion" : 0e+10,
                                                "rotation_2" : 0e+10,
                                                "rotation_3" : 0e+10})

# from KratosMultiphysics import *
# from KratosMultiphysics.IgaApplication import *
# # beam_a.nodes[0].Fix(DISPLACEMENT_X)
# # beam_a.nodes[0].Fix(DISPLACEMENT_Y)
# # beam_a.nodes[0].Fix(DISPLACEMENT_Z)
# # beam_a.nodes[0].Fix(DISPLACEMENT_ROTATION)
# beam_a.nodes[-1].Fix(DISPLACEMENT_X)
# beam_a.nodes[-1].Fix(DISPLACEMENT_Y)
# beam_a.nodes[-1].Fix(DISPLACEMENT_Z)
# beam_a.nodes[-1].Fix(DISPLACEMENT_ROTATION)

beam_a.add_node_load(index=-1)
# beam_a.add_node_load_moment(t=beam_a.t1(), force=[0, 0, 0], moment=[0,0,0])
# beam_b.add_node_load(index=-1)

model.init_solver()

for step, lam in enumerate(np.linspace(0, 1, 6)[1:]):
    beam_a.make_header(step)

    print(' ##############################')
    print('Force lambda-z: ', lam )
    F = 2 * lam
    print(F)

    # beam_a.set_node_value(index=-1, directions=['rotation'], value=F)

    # beam_a.set_node_load(index=-1, load=[ -F ,0, 0])
    beam_a.set_node_load(index=-1, load=[0, -F, 0])
    # w = (F * 4**3)/(48 * youngs_modulus * 1)     # Einfeldträger
    # w = (F*10)**2/(youngs_modulus * 1)          # Kragarmpaper
    # print('w = ' + '%.12f' % (w) )      

    # beam_b.set_node_load(index=-1, load=[0,0 , -0.00001*lam ])

    model.solve(lam)

    beam_a.write_displacement(step)
    beam_b.write_displacement(step)
    beam_a.write_results(step)
    geometry.Save(r'data/output.iga')


# --- output

beam_a.print_forces(1)
print('BEAM_A')
beam_a.print_displacement()

geometry.Save(r'data/output.iga')

print('Done!')
