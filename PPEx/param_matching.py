"""Based on the laser and plasma parameters, computes relevant physical quantities."""


from math import pi, log, sqrt


class Laser(object):
    """bi-Gaussian laser beam propagating along the z axis."""
    
    c = 0.2998 # speed of light in vaccum in microns/fs

    def __init__(self, w_0, lambda_0, tau_0, P=None, epsilon=None):
        """Create laser with given parameters.

        w_0 -- beam waist in microns at the focal plane z = 0
        lambda_0 -- wavelength in microns
        tau_0 -- FWHM of the pulse duration in fs

        One can give either
        P -- laser power in TW --or--
        epsilon -- laser energy in Joules
        """
        self.w_0 = w_0
        self.lambda_0 = lambda_0
        self.tau_0 = tau_0
        self.k_0 = 2 * pi / lambda_0 # wavenumber in rads/micron
        self.omega_0 = Laser.c * self.k_0  # angular frequency in rads/fs
        self.z_R = pi * self.w_0**2 / self.lambda_0 # Rayleigh length in microns

        if epsilon and not P:
            self.epsilon = epsilon
            self.P = self.power()
        elif P and not epsilon:
            self.P = P
            self.epsilon = self.energy()
        else:
            raise TypeError('Either laser power or energy must be given.')
        
        # peak intensity in the focal plane, in 10^20 W/cm^2
        self.I_0 = 2 * self.P / (pi * self.w_0**2)
        # rescaled vector potential, dimensionless
        self.a_0 = 8.6 * self.lambda_0 * sqrt(self.I_0)
        # amplitude of the electric field in TV/m
        self.E = 3.21 * self.a_0 / self.lambda_0


    def power(self):
        """Compute power, given pulse energy and duration."""
        pw_to_tw = 1e3
        return 2 * sqrt(log(2) / pi) * self.epsilon / self.tau_0 * pw_to_tw
    
    def energy(self):
        """Given power, compute the pulse energy."""
        tw_to_pw = 1e-3
        return self.P * tw_to_pw * self.tau_0 / (2 * sqrt(log(2) / pi))



class Plasma(object):
    """Electron plasma with fixed density."""

    def __init__(self, n_p):
        """Create plasma with given density.
        
        n_e -- plasma density in units of 10^18 cm^(-3)
        """
        self.n_p = n_p
        self.lambda_p = 33.39 / sqrt(self.n_p) # skin depth in microns
        self.k_p = 2 * pi / self.lambda_p # wavenumber in rad/micron
        self.omega_p = 5.64 * 1e-2 * sqrt(self.n_p) # electron plasma frequency in rad/fs


class Matching(object):
    """Matching conditions for laser-plasma interaction"""

    def __init__(self, laser, plasma):
        """Generate conditions based on certain laser and plasma parameters
        
        laser -- object of the Laser class
        plasma -- object of the Plasma class
        """
        self.laser = laser
        self.plasma = plasma
        # critical density in units of 10^18 cm^(-3)
        self.n_c = 1120 / self.laser.lambda_0**2
        
        omega_ratio = (self.laser.omega_0 / self.plasma.omega_p)**2
        # dephasing length in microns
        self.L_d = 4 / 3 * omega_ratio * sqrt(self.laser.a_0) / self.plasma.k_p
        # pump depletion length in microns
        self.L_pd = omega_ratio * Laser.c * self.laser.tau_0
        # critical power in TW for relativistic self-focusing
        self.P_c = 0.017 * omega_ratio

    def conditions(self):
        message = "Current plasma density is n_p = {:06.2f} 10^18 cm^(-3).\n".format(self.plasma.n_p)
        waist = sqrt(self.laser.a_0) * self.plasma.lambda_p / pi
        message += "The recommended beam waist at this density is {:06.2f} micron.\n".format(waist)
        message += "The recommended a_0 = {:06.2f}.\n".format(2* (self.laser.P / self.P_c)**(1/3))
        if self.laser.a_0 >= 4:
            nlim = self.n_c / self.laser.a_0**5
            message += "Etching rate should exceed diffraction rate: n_p >= {:06.2f} 10^18 cm^(-3).\n".format(nlim)
        nlim = 30 / self.laser.P
        message += "For good self-guiding, P >= P_c, therefore n_p >= {:06.2f} 10^18 cm^(-3).\n".format(nlim)
        nlim = 112.96 * self.laser.a_0 / self.laser.w_0**2
        message += "To match the bubble radius to the beam waist, n_p <= {:06.2f} 10^18 cm^(-3).\n".format(nlim)
        return message



if __name__ == "__main__":
    # Lu parameters
    # laser = Laser(w_0=19.5, lambda_0=0.8, tau_0=30, P=200)
    # plasma = Plasma(n_p=1.5)

    # CETAL parameters
    # laser = Laser(w_0=30, lambda_0=0.8, tau_0=40, epsilon=7)
    # plasma = Plasma(n_p=0.294)

    # CALDER-CIRC_2008.pdf
    laser = Laser(w_0=18, lambda_0=0.8, tau_0=30, P=17.303)
    plasma = Plasma(n_p=7.5)
    
    # FBPIC lwfa_script.py
    # laser = Laser(w_0=5, lambda_0=0.8, tau_0=16.67, P=13.27)
    # plasma = Plasma(n_p=4)


    print("Laser parameters:")
    print('Rayleigh length z_R = {:04.2f} mm.'.format(laser.z_R*1e-3))
    print('Peak intensity in the focal plane: I_0 = {:03.3f} 10^20 W/cm^2'.format(laser.I_0))
    print('a_0 = {:04.2f}'.format(laser.a_0))
    print('Electric field E_L = {:04.2f} TV/m.'.format(laser.E))
    print()


    
    print('Plasma parameters:')
    print('Plasma skin depth lambda_p = {:04.2f} micron'.format(plasma.lambda_p))
    print('Plasma wave number k_p = {:04.2f} rad/micron'.format(plasma.k_p))
    print('Plasma frequency omega_p = {:04.3f} rad/fs'.format(plasma.omega_p))
    print()

    m = Matching(laser, plasma)
    print('Matching conditions:')
    print('Dephasing length L_d = {:08.1f} mm.'.format(m.L_d*1e-3))
    print('Pump depletion length L_pd = {:08.1f} mm.'.format(m.L_pd*1e-3))
    print('Critical power for relativistic self-focusing P_c = {:04.2f} TW.'.format(m.P_c))
    print()
    print(m.conditions())