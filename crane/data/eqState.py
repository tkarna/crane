#!/usr/bin/python
"""
Implementation of Jackett et al. (2006) equation of state.

Vertical density gradient is best computed with the thermal expansion and haline
contraction coefficients (alpha and beta, respectively):
dRhodz = alpha(S,Th,p)*dThdz_at_const_p + beta(S,Th,p)*dSdz_at_const_p

Jackett, D. R., McDougall, T. J., Feistel, R., Wright, D. G., and Griffies, S. M. (2006).
Algorithms for Density, Potential Temperature, Conservative Temperature, and the Freezing Temperature of Seawater. Journal of Atmospheric and Oceanic Technology, 23(12):1709-1728.

Tuomas Karna 2013-07-31
"""

import numpy as np

# Coefficients for the rational function from Jackett et al. (2006)
a = [9.9984085444849347e2, 7.3471625860981584e0, -5.3211231792841769e-2,
     3.6492439109814549e-4, 2.5880571023991390e0, -6.7168282786692355e-3,
     1.9203202055760151e-3, 1.1798263740430364e-2, 9.8920219266399117e-8,
     4.6996642771754730e-6, -2.5862187075154352e-8, -3.2921414007960662e-12]
b = [1.0, 7.2815210113327091e-3, -4.4787265461983921e-5, 3.3851002965802430e-7,
     1.3651202389758572e-10, 1.7632126669040377e-3, -8.8066583251206474e-6,
     -1.8832689434804897e-10, 5.7463776745432097e-6, 1.4716275472242334e-9,
     6.7103246285651894e-6, -2.4461698007024582e-17, -9.1534417604289062e-18]


def nominator(S, Th, p):
    """Computes the nominator of the equation of state."""
    Pn = (
        a[0] +
        Th *
        a[1] +
        Th *
        Th *
        a[2] +
        Th *
        Th *
        Th *
        a[3] +
        S *
        a[4] +
        Th *
        S *
        a[5] +
        S *
        S *
        a[6] +
        p *
        a[7] +
        p *
        Th *
        Th *
        a[8] +
        p *
        S *
        a[9] +
        p *
        p *
        a[10] +
        p *
        p *
        Th *
        Th *
        a[11])
    return Pn


def denominator(S, Th, p):
    """Computes the denominator of the equation of state."""
    Pd = (b[0] +
          Th *
          b[1] +
          Th *
          Th *
          b[2] +
          Th *
          Th *
          Th *
          b[3] +
          Th *
          Th *
          Th *
          Th *
          b[4] +
          S *
          b[5] +
          S *
          Th *
          b[6] +
          S *
          Th *
          Th *
          Th *
          b[7] +
          pow(S, 1.5) *
          b[8] +
          pow(S, 1.5) *
          Th *
          Th *
          b[9] +
          p *
          b[10] +
          p *
          p *
          Th *
          Th *
          Th *
          b[11] +
          p *
          p *
          p *
          Th *
          b[12])
    return Pd


def equationOfStateJackett(S, Th, p):
    """Computes seawater density, based on salinity (S), potential temperature (Th) and pressure (p)."""
    Pn = nominator(S, Th, p)
    Pd = denominator(S, Th, p)
    rho = Pn / Pd
    return rho


def equationOfStateJackettBeta(S, Th, p):
    """Computes seawater haline contraction coeffcient beta, based on salinity (S), potential temperature (Th) and pressure (p)."""
    Pn = nominator(S, Th, p)
    dndS = a[4] + a[5] * Th + 2 * a[6] * S + a[9] * p
    Pd = denominator(S, Th, p)
    dddS = (b[5] + b[6] * Th + b[7] * Th * Th * Th + 1.5 * b[8] * pow(S, 0.5) +
            1.5 * b[9] * pow(S, 0.5) * Th * Th)
    beta = (Pd * dndS - Pn * dddS) / Pd / Pd
    return beta


def equationOfStateJackettAlpha(S, Th, p):
    """Computes seawater thermal expansion coeffcient alpha, based on salinity (S), potential temperature (Th) and pressure (p)."""
    Pn = nominator(S, Th, p)
    dndT = a[1] + 2 * a[2] * Th + 3 * a[3] * Th * Th + a[
        5] * S + 2 * a[8] * p * Th + 2 * a[11] * p * p * Th
    Pd = denominator(S, Th, p)
    dddT = (b[1] +
            2 *
            b[2] *
            Th +
            3 *
            b[3] *
            Th *
            Th +
            4 *
            b[4] *
            Th *
            Th *
            Th +
            b[6] *
            S +
            3 *
            b[7] *
            S *
            Th *
            Th +
            2 *
            b[9] *
            pow(S, 1.5) *
            Th +
            3 *
            b[11] *
            p *
            p *
            Th *
            Th +
            b[12] *
            p *
            p *
            p)
    alpha = (Pd * dndT - Pn * dddT) / Pd / Pd
    return alpha
