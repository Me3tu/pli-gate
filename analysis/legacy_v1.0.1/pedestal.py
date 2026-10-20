#!/usr/bin/env python3
"""Truncated-exponential mean, its inversion, and the over-subtraction bias.

A delayed-coincidence lifetime measurement keeps events in a finite window
(tau_min, T). The mean of a truncated exponential is not the lifetime, so the
observed mean has to be inverted to recover tau. If the randoms estimator
removes more counts than it should, the residue behaves as a uniform pedestal
across the window and biases the recovered lifetime.

    m1(tau)      mean of Exp(tau) truncated to (tau_min, T), minus tau_min
    invert_trunc inverse of m1: recovers tau from an observed mean
    predicted_bias(tau, phi)
                 fractional error in tau-hat caused by over-subtracting a
                 fraction phi of uniform counts:

                     tau_obs(phi) = [m1(tau) - phi*m_ped] / (1 - phi)
                     m_ped        = (T - tau_min)/2          (uniform mean)
                     bias         = tau_hat(phi)/tau_hat(0) - 1

phi is a measured quantity, not a fitted one.
"""
from __future__ import annotations
import math

TAU_MIN = 1.0     # ns, lower censor
T_HI = 5.0        # ns, upper edge of the prompt window
WIDTH = T_HI - TAU_MIN


def trunc_mean(tau: float, lo: float = TAU_MIN, hi: float = T_HI) -> float:
    """Mean of Exp(tau) truncated to (lo, hi), measured from lo."""
    if tau <= 0:
        return float("nan")
    a, b = lo / tau, hi / tau
    ea, eb = math.exp(-a), math.exp(-b)
    num = (lo + tau) * ea - (hi + tau) * eb
    den = ea - eb
    return num / den - lo


def invert_trunc(m_target: float, lo: float = TAU_MIN, hi: float = T_HI,
                 tol: float = 1e-12) -> float:
    """Recover tau from a truncated mean by bisection. Monotone, so this is safe."""
    if not (0 < m_target < hi - lo):
        return float("nan")
    a, b = 1e-6, 1e4
    for _ in range(400):
        mid = 0.5 * (a + b)
        if trunc_mean(mid, lo, hi) < m_target:
            a = mid
        else:
            b = mid
        if b - a < tol:
            break
    return 0.5 * (a + b)


def predicted_bias(tau_true: float, phi: float) -> float:
    """Fractional bias in tau-hat from over-subtracting phi*N uniform counts."""
    if phi >= 1.0:
        return float("nan")
    m1 = trunc_mean(tau_true)
    m_ped = WIDTH / 2.0
    m_obs = (m1 - phi * m_ped) / (1.0 - phi)
    if not (0 < m_obs < WIDTH):
        return float("nan")
    return invert_trunc(m_obs) / invert_trunc(m1) - 1.0
