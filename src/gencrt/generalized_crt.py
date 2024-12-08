"""
generalized_crt.py

This module implements the generalized Chinese remainder theorem (CRT)
algorithm.

Author: Simon Ljungbeck
Date: December 2024
"""

from typing import Iterable, Tuple


def extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """
    Compute the greatest common divisor (GCD) of two integers a and b using
    the Extended Euclidean Algorithm. It also computes the coefficients of
    Bézout's identity, i.e., integers x and y such that:
    a * x + b * y = gcd(a, b).

    Args:
        a (int): The first integer.
        b (int): The second integer.

    Returns:
        Tuple[int, int, int]: A tuple containing:
            - gcd (int): The greatest common divisor of a and b.
            - x (int): The coefficient of a in Bézout's identity.
            - y (int): The coefficient of b in Bézout's identity.
    """
    if b == 0:
        return a, 1, 0
    gcd, x1, y1 = extended_gcd(b, a % b)
    x = y1
    y = x1 - (a // b) * y1
    return gcd, x, y


def crt_two_eqs(a: int, m: int, b: int, n: int) -> int:
    """
    Compute the solution to a modular equation system with two equations:
    x ≡ a (mod m) and x ≡ b (mod n).

    Args:
        a (int): The remainder for the first equation.
        m (int): The modulus for the first equation.
        b (int): The remainder for the second equation.
        n (int): The modulus for the second equation.

    Returns:
        int: The solution x modulo lcm(m, n) to the equation system, or None
        if a solution does not exist.
    """
    g, u, v = extended_gcd(m, n)
    if a % g == b % g:
        return ((a * v * n + b * u * m) // g) % (m * n)
    return None


def crt(congruences: Iterable[Tuple[int, int]]) -> int:
    """
    Compute the solution to a modular equation system using the generalized
    Chinese remainder theorem (CRT).

    Args:
        congruences (Iterable[Tuple[int, int]]): A list of tuples (a_i, n_i)
        such that x ≡ a_i (mod n_i) for all (a_i, n_i).

    Returns:
        int: The solution x modulo prod(n_i) to the equation system, or None
        if a solution does not exist.
    """
    if len(congruences) == 0:
        return None
    if len(congruences) == 1:
        return congruences[0][0] % congruences[0][1]
    a0, n0 = congruences[0]
    for a1, n1 in congruences[1:]:
        x = crt_two_eqs(a0, n0, a1, n1)
        if x is None:
            return None
        a0 = x
        n0 *= n1
    return a0
