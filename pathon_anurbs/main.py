from kratos_utilities import Model
import os
import ANurbsDev as an
import numpy as np

print('Process ID: ', os.getpid())
print(' ')

geometry = an.Model()
geometry.Load(r'data/model_curve.iga')

model = Model(geometry)

for curve_geometry_ptr in geometry.GetByType('CurveGeometry3D'):
    model.add_beam(curve_geometry_ptr)

a = 0.01       # Querschnittshöhe
b = 0.01       # Querschnittsbreite
c1 = 1/3 * ( 1 - 0.63/(a/b) + 0.052/((a/b)**5))
youngs_modulus = 1.999e8

model.add_beam_properties('material_1',
    area = a * b,   # Querschnittsfläche
    it = c1 * a * b**3,
    iy = (a * b**3)/12,   # Flächenträgheitsmoment Iy
    iz = (b * a**3)/12,   # Flächenträgheitsmoment Iz
    youngs_modulus = youngs_modulus,
    shear_modulus = youngs_modulus / 2,
)


beam_a = model.beam('curve_a')
# beam_b = model.beam('curve_b')

beam_a.add_stiffness(
    material='material_1',
)

# beam_b.add_stiffness(
#     material='material_1',
# )

# beam_a.add_support(t=beam_a.t0)
# beam_a.add_coupling(t=beam_a.t1, other=beam_b, other_t=beam_b.t0)
# beam_b.add_support(t=beam_a.t1)


# beam_a.fix_node(index=0, directions=['x', 'y', 'z', 'rotation'])
# beam_a.fix_node(index=1, directions=['y', 'z'])

penalty_dict = {"displacement" : 1e+7, 
                "torsion" : 1e+7,
                "rotation" : 1e7}

beam_a.add_support(t = 0, penalty = penalty_dict, material='material_1')

beam_a.add_node_load(index=-1)



for step, lam in enumerate(np.linspace(0, 1, 10)[1:]):
    beam_a.set_node_load(index=-1, load=[0, 0, 0.001*lam])
    model.solve()

# --- output

geometry.Save(r'data/output.iga')

print('Done!')