from kratos_utilities import Model
import os
import ANurbsDev as an
import numpy as np

print('Process ID: ', os.getpid())
print(' ')

geometry = an.Model()
geometry.Load(r'data/einfeldtraeger.iga')
open('OutputAD.txt', 'w').close()
open('normalkraft.txt', 'w').close()
open('moment_2.txt', 'w').close()
open('moment_3.txt', 'w').close()


model = Model(geometry)

for curve_geometry_ptr in geometry.GetByType('CurveGeometry3D'):
    model.add_beam(curve_geometry_ptr)

a = 0.01       # Querschnittshöhe
b = 0.01       # Querschnittsbreite
c1 = 1/3 * ( 1 - 0.63/(a/b) + 0.052/((a/b)**5))
youngs_modulus = 100

model.add_beam_properties('material_1',
    # area = a * b,   # Querschnittsfläche
    # it = c1 * a * b**3,
    # iy = (a * b**3)/12,   # Flächenträgheitsmoment Iy
    # iz = (b * a**3)/12,   # Flächenträgheitsmoment Iz
    area = 20000000,   # Querschnittsfläche
    it = 0.5,
    iy = 1,   # Flächenträgheitsmoment Iy
    iz = 1,   # Flächenträgheitsmoment Iz
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


# beam_a.fix_node(index=-1, directions=['x', 'y', 'z'])
# beam_a.fix_node(index=1, directions=['y', 'z'])

beam_a.add_support(t = beam_a.t0(), penalty = {"displacement_x" : 1e+10,
                                                "displacement_y" : 1e+10,
                                                "displacement_z" : 1e+10, 
                                                "torsion" : 0e+10,
                                                "rotation_2" : 0e+10,
                                                "rotation_3" : 0e+10})

beam_b.add_support(t = beam_b.t1(), penalty = {"displacement_x" : 1e+10,
                                                "displacement_y" : 0e+10,
                                                "displacement_z" : 1e+10, 
                                                "torsion" : 0e+10,
                                                "rotation_2" : 0e+10,
                                                "rotation_3" : 0e+10})


beam_a.write_gausspoints()


beam_a.add_node_load(index=-1)
# beam_b.add_node_load(index= -1)

for step, lam in enumerate(np.linspace(0, 1, 11)[1:]):
    print(' ##############################')
    print('Force lambda-z: ', lam )
    F = 15*lam
    # beam_a.set_node_load(index=-1, load=[-0.001*lam, 0 , 0])
    # beam_a.set_node_load(index=-1, load=[0,  50*lam ,0])
    beam_a.set_node_load(index=-1, load=[0, 0 , F])
    w = (F * 4**3)/(48 * youngs_modulus * 1)     # Einfeldträger
    # w = (F*10)**2/(youngs_modulus * 1)          # Kragarmpaper
    print('w = ' + '%.12f' % (w) )      

    # beam_b.set_node_load(index=-1, load=[0,0 , 1*lam ])

    # moment = [0, 1*lam, 0]

    # beam_a.add_moment(t = beam_a.t1(), vector = moment, material = 'material_1' )

    model.solve()

# --- output

beam_a.print_forces()

geometry.Save(r'data/output.iga')

print('Done!')

beam_a.print_displacement()