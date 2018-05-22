import pytest
import param_matching as pm
from pytest import approx

def test_laser():
    laser = pm.Laser(w_0=30, lambda_0=0.8, tau_0=40, epsilon=7)
    assert laser.z_R == approx(3530, 1)
    assert laser.k_0 == approx(7.853, 1e-3)
    assert laser.omega_0 == approx(2.3546, 1e-3)
    assert laser.P == approx(164.4, 0.1)
    assert pm.Laser(w_0=30, lambda_0=0.8, tau_0=40, P=164.4).epsilon == approx(7, 1e-5)
    with pytest.raises(TypeError):
        pm.Laser(w_0=30, lambda_0=0.8, tau_0=40)
    assert laser.I_0 == approx(0.116, 1e-2)
    assert laser.a_0 == approx(2.343, 1e-2)
    assert laser.E == approx(9.4, 0.1)

def test_plasma():
    pl = pm.Plasma(n_p=0.294)
    assert pl.lambda_p == approx(61.58, 1e-2)
    assert pl.k_p == approx(0.102, 1e-2)
    assert pl.omega_p == approx(0.0305, 1e-2)

def test_matching():
    laser = pm.Laser(w_0=30, lambda_0=0.8, tau_0=40, epsilon=7)
    plasma = pm.Plasma(n_p=0.294)
    m = pm.Matching(laser, plasma)
    assert m.n_c == approx(1750, 1)
    assert m.L_d == approx(118664, 1)
    assert m.L_pd == approx(71093, 1)
    assert m.P_c == approx(100, 1)    

