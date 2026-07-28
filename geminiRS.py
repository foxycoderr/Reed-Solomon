# Change this:
# from numpy import polynomial

# To this:
from numpy.polynomial.polynomial import Polynomial as polynomial

import numpy as isaac_newton
import random
from sympy import Matrix

alpha = 6
mod = 251


def polymul_mod(poly1, poly2, modulus=251):
    """Multiplies two polynomials and keeps coefficients modulo 251 to avoid float bugs."""
    c1 = [int(round(c)) % modulus for c in poly1.coef]
    c2 = [int(round(c)) % modulus for c in poly2.coef]

    res = [0] * (len(c1) + len(c2) - 1)
    for i, a in enumerate(c1):
        for j, b in enumerate(c2):
            res[i + j] = (res[i + j] + a * b) % modulus

    return polynomial.Polynomial(res)


def encode(n, k, message, alpha=6):
    if len(message) != k:
        print(f"Invalid message length ({len(message)}).")
        exit(1)

    for message_symbol in message:
        if not isinstance(message_symbol, int) or message_symbol > 250 or message_symbol < 0:
            print("Your message is cooked")
            exit(1)

    message_polynomial = polynomial.Polynomial(message)
    generator_polynomial = polynomial.Polynomial([-int(alpha) % 251, 1])

    for i in range(2, n - k + 1):
        next_root_poly = polynomial.Polynomial([-pow(int(alpha), i, 251), 1])
        generator_polynomial = polymul_mod(generator_polynomial, next_root_poly, 251)

    codeword_polynomial = polymul_mod(generator_polynomial, message_polynomial, 251)

    # Pad with trailing zeros up to block size n
    coeffs = [int(round(c)) % 251 for c in codeword_polynomial.coef]
    while len(coeffs) < n:
        coeffs.append(0)

    return polynomial.Polynomial(coeffs)


def transmit(codeword_poly, e, n=10):
    """Safely injects exactly e errors without ruining array structure."""
    coeffs = [int(round(c)) % 251 for c in codeword_poly.coef]
    while len(coeffs) < n:
        coeffs.append(0)

    for _ in range(e):
        index_of_error_location = random.randint(0, len(coeffs) - 1)
        error_magnitude = random.randint(1, 250)
        coeffs[index_of_error_location] = (coeffs[index_of_error_location] + error_magnitude) % 251

    return polynomial.Polynomial(coeffs)


def find_syndromes(n, k, r_polynomial):
    """Evaluates syndromes using pure modular Horner's method to avoid huge floats."""
    syndromes = [0]  # 0-indexed padding so S1 matches index 1
    coeffs = [int(round(c)) % 251 for c in r_polynomial.coef]

    for j in range(1, n - k + 1):
        Sj = 0
        root = pow(alpha, j, 251)
        # Horner's scheme for modular evaluation
        for coeff in reversed(coeffs):
            Sj = (Sj * root + coeff) % 251
        syndromes.append(Sj)
    return syndromes


def find_elp(syndromes, p_max):
    """Finds ELP dynamically by finding the actual error count up to p_max."""
    p = p_max
    while p > 0:
        syndromes_pxp_matrix_raw = []
        for i in range(1, p + 1):
            row = []
            for j in range(0, p):
                row.append(syndromes[i + j])
            syndromes_pxp_matrix_raw.append(row)

        syndomes_pxp_matrix_sympy = Matrix(syndromes_pxp_matrix_raw)

        # If determinant is zero, fewer errors happened. Step down!
        if syndomes_pxp_matrix_sympy.det() % 251 == 0:
            p -= 1
            continue

        inverse_syndromes_pxp_matrix_sympy = syndomes_pxp_matrix_sympy.inv_mod(251)
        inverse_syndromes_pxp_matrix = isaac_newton.matrix(inverse_syndromes_pxp_matrix_sympy)

        syndomes_px1_matrix_raw = []
        for i in range(p + 1, 2 * p + 1):
            syndomes_px1_matrix_raw.append([-syndromes[i]])
        syndromes_px1_matrix = isaac_newton.matrix(syndomes_px1_matrix_raw)

        elp_coefficients = isaac_newton.mod(inverse_syndromes_pxp_matrix * syndromes_px1_matrix, 251)
        elp_coefficients = [int(c[0, 0]) for c in list(elp_coefficients)]
        elp_coefficients.insert(0, 1)

        return polynomial.Polynomial(elp_coefficients), p

    return polynomial.Polynomial([1]), 0  # 0 Errors found


def solve_elp(error_correcting_polynomial):
    Xr_inverses = []
    # Test all 251 elements of the field
    for i in range(251):
        raw_val = error_correcting_polynomial(i)
        if abs(round(raw_val) - raw_val) < 1e-7:
            if round(raw_val) % 251 == 0:
                Xr_inverses.append(i)

    # Convert inverses to error locations Xr = 1 / Xr_inverse
    fermat_inverse_vec = isaac_newton.vectorize(lambda a, p: pow(int(a), int(p - 2), int(p)))
    if not Xr_inverses:
        return []
    Xr = list(fermat_inverse_vec(Xr_inverses, 251))
    return [int(i) for i in Xr]


def find_error_polynomial_coefficients_Yr(syndromes, Xr, p, n, k):
    if p == 0 or not Xr:
        return []

    augmented_matrix = []
    # Constructing Syndromes parity matrix rows for the error locations
    for row_power in range(1, n - k + 1):
        row = []
        for X in Xr:
            row.append(pow(X, row_power, 251))
        row.append(syndromes[row_power])
        augmented_matrix.append(row)

    # Clean Row Reduction (Gaussian elimination) mod 251
    num_vars = len(Xr)
    for i in range(num_vars):
        pivot = augmented_matrix[i][i]
        pivot_inverse = pow(pivot, -1, 251)

        for j in range(num_vars + 1):
            augmented_matrix[i][j] = (augmented_matrix[i][j] * pivot_inverse) % 251

        for target_row in range(num_vars):
            if i == target_row:
                continue
            factor = augmented_matrix[target_row][i]
            for j in range(num_vars + 1):
                augmented_matrix[target_row][j] = (augmented_matrix[target_row][j] - factor * augmented_matrix[i][
                    j]) % 251

    Yr = []
    for i in range(num_vars):
        Yr.append(augmented_matrix[i][num_vars])
    return Yr


def generate_error_polynomial_coefficients(Xr, Yr, n):
    error_polynomial_coefficients = [0] * n
    for X, Y in zip(Xr, Yr):
        # Locate the exact index i where alpha^i == X
        for i in range(251):
            if pow(alpha, i, 251) == X:
                error_polynomial_coefficients[i] = Y
                break
    return error_polynomial_coefficients


def decode(n, k, received):
    r_polynomial = polynomial.Polynomial(received)
    syndromes = find_syndromes(n, k, r_polynomial)
    t = int(isaac_newton.floor((n - k) / 2))

    error_correcting_polynomial, actual_p = find_elp(syndromes, t)

    if actual_p == 0:
        # Zero errors found: clean message directly
        return [int(round(c)) % 251 for c in r_polynomial.coef][:k]

    Xr = solve_elp(error_correcting_polynomial)
    Yr = find_error_polynomial_coefficients_Yr(syndromes, Xr, actual_p, n, k)

    error_coeffs = generate_error_polynomial_coefficients(Xr, Yr, n)
    error_polynomial = polynomial.Polynomial(error_coeffs)

    message_polynomial = polymul_mod(r_polynomial - error_polynomial, polynomial.Polynomial([1]), 251)

    # Return message coefficients stripped back to k length
    clean_coeffs = [int(round(c)) % 251 for c in message_polynomial.coef]
    return clean_coeffs[:k]


# Execution Test
codeword = encode(10, 5, [1, 2, 3, 4, 5])
print("Encoded Codeword:", [int(c) for c in codeword.coef])

# Test 1: Direct decoding with ZERO errors
decoded_clean = decode(10, 5, codeword)
print("Decoded (0 errors):", decoded_clean)

# Test 2: Injecting 1 random error
corrupted_codeword = transmit(codeword, e=1)
decoded_corrupted = decode(10, 5, corrupted_codeword)
print("Decoded (1 error): ", decoded_corrupted)
