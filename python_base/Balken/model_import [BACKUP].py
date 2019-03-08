import ANurbs as an

model = an.Model.open(r'C:\_Masterarbeit\beispiele\Balken\Balken.iga')

curve_item = model.of_type('Curve3D')[0]
curve = curve_item.data
curve_geometry = curve_item.geometry().geometry

print('Kontrollpunkte und Gewichte:')

for i in range(curve_geometry.NbPoles):
    pole = curve_geometry.Pole(i)
    weight = curve_geometry.Weight(i)

    print('  ', i, pole, weight)

print()
print('Knots:')
knots = curve_geometry.Knots
print('  ', knots)
print()
print('Integrationspunkte und Gewichte im Parameterraum der Kurve:')

integration_points = curve_item.IntegrationPoints()

for t, weight in integration_points:
    print('  ', 't =', t, 'weight =', weight)

print()
print('Tangentenvektor:')

integration_points = curve_item.IntegrationPoints()

for t, weight in integration_points:
    point, tangent = curve.DerivativesAt(T=t, Order=1)
    print('  ', 't =', t, 'tangent =', tangent)

print()
print('Formfunktionen:')

shapes = an.CurveShapeEvaluator(Degree=curve_geometry.Degree, Order=2)

for t, weight in integration_points:
    shapes.Compute(curve_geometry.Knots, t)

    n_0 = [0] * shapes.NbNonzeroPoles
    n_1 = [0] * shapes.NbNonzeroPoles
    n_2 = [0] * shapes.NbNonzeroPoles

    for i in range(shapes.NbNonzeroPoles):
        n_0[i] = shapes(0, i)
        n_1[i] = shapes(1, i)
        n_2[i] = shapes(2, i)

    print('  ', 't =', t)
    print('    ', 'n_0 =', n_0)
    print('    ', 'n_1 =', n_1)
    print('    ', 'n_2 =', n_2)

print(curve_geometry.PointAt(1.0))
print('done!')