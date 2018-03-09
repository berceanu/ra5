import deck_classes as dc
import pint
import math

ureg = pint.UnitRegistry()
ureg.load_definitions('./units/lwfa_def.txt')

mmovingwin = dc.MovingWindow(dc.WindowProp(0.29919287308 * ureg.micrometer/ureg.femtosecond, 166.78204759999997 * ureg.femtosecond, 9041.7551465 * ureg.femtosecond),
                          dc.Point(-100 * ureg.micrometer, -400 * ureg.micrometer, -400 * ureg.micrometer),
                          dc.Point(-50 * ureg.micrometer, 400 * ureg.micrometer, 400 * ureg.micrometer),
                          dc.NrGridPoints(1024, 96, 96))
ddomain = dc.Domain(mmovingwin)
pplasma = dc.Plasma(2.6575 * ureg.micrometer, 100. * ureg.micrometer, 2450. * ureg.micrometer)
eelectron = dc.Species(ddomain, pplasma, 1e-3)

def test_domain():
    assert ddomain.mw.width.x == 50.0 * ureg.micrometer
    assert ddomain.mw.width.y / 2.0 == 400.0 * ureg.micrometer
    assert ddomain.mw.box.x.shape == (1024, )
    assert ddomain.npoints.x == 55352 and math.isclose(ddomain.width.x.magnitude, 2705.3286999676016)

def test_plasma():
    assert math.isclose(pplasma.density.magnitude, 3.9986353164000005e+18, rel_tol=1e-5)

def test_rho_min():
    assert math.isclose(eelectron.rho_min.magnitude, 3998635316400000.5, rel_tol=1e-5)
