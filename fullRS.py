from numpy import polynomial
import numpy as isaac_newton
import random

from sympy import Matrix

# Globals
alpha = 6
mod = 251

# ENCODER -- Debugged.
def polymul_mod(poly1, poly2, modulus=251):
    """Multiplies two polynomials and keeps coefficients modulo 251 to avoid float bugs."""
    # Convert inputs to clean lists of raw integers
    c1 = [int(round(c)) % modulus for c in poly1.coef]
    c2 = [int(round(c)) % modulus for c in poly2.coef]

    # Perform manual integer-safe multiplication (convolution)
    res = [0] * (len(c1) + len(c2) - 1)
    for i, a in enumerate(c1):
        for j, b in enumerate(c2):
            res[i + j] = (res[i + j] + a * b) % modulus

    return polynomial.Polynomial(res)
def encode(n, k, message):
    if len(message) != k:
        print(f"Invalid message length ({len(message)}).")
        exit(1)

    for message_symbol in message:
        if not isinstance(message_symbol, int) or message_symbol > 250 or message_symbol < 0:
            print("Your message is cooked")
            exit(1)

    # 1. Initiate message and the first root of the generator polynomial
    message_polynomial = polynomial.Polynomial(message)
    generator_polynomial = polynomial.Polynomial([-int(alpha) % 251, 1])

    # 2. Safely chain multiplications using our modular helper function
    for i in range(2, n - k + 1):
        next_root_poly = polynomial.Polynomial([-pow(alpha, i, 251), 1])
        generator_polynomial = polymul_mod(generator_polynomial, next_root_poly, 251)

    # 3. Create final codeword safely modulo 251
    codeword_polynomial = polymul_mod(generator_polynomial, message_polynomial, 251)

    print("Clean Codeword Polynomial:\n", codeword_polynomial)
    return codeword_polynomial

# ERROR INTRODUCTION - Debugged.
def transmit(message, e):
    message = list(message)
    error_locations = []
    error_values = []
    for i in range(e):
        index_of_error_location = random.randint(0, len(message) - 1)
        error_locations.append(index_of_error_location)
        error_values.append(random.randint(1, 250))
        message[index_of_error_location] += error_values[-1]
        message[index_of_error_location] = index_of_error_location % 251

    message = polynomial.Polynomial(message)

    return [message, error_locations, error_values]

# DECODER
def eval_mod_poly(coeffs, x, p=251):
    y = 0
    for c in reversed(coeffs):
        y = (y * x + c) % p
    return y

def find_syndromes(n, k, r_polynomial):
    print("R-poly clean int coeffs in find_syndromes")
    r_coeffs = list(int(i) for i in list(r_polynomial.coef[0]))
    print(r_coeffs)
    print()

    syndromes = []
    for j in range(1, n - k + 1):
        x = pow(alpha, j, mod)
        Sj = eval_mod_poly(r_coeffs, x, mod)
        print(f"S{j}", Sj)
        syndromes.append(Sj)
        print()
    print("Syndromes", syndromes)
    return syndromes

def find_elp(syndromes, t):
    possible_error_count = reversed(list(i for i in range(t+1)))
    inverse_syndromes_pxp_matrix = None
    print("Beginning search for pxp matrix. ")
    for p in possible_error_count:
        # finding large LHS matrix
        print(">>> Trying p = ", p, "...")
        print()
        syndromes_pxp_matrix_raw = []
        for i in range(1, p+1):
            row = []
            for j in range(0, p):
                row.append(syndromes[i+j - 1])
            syndromes_pxp_matrix_raw.append(row)

        print(f"Syndromes {p}x{p} matrix raw:")
        print(syndromes_pxp_matrix_raw)
        # inverting LHS matrix
        syndomes_pxp_matrix = isaac_newton.matrix(syndromes_pxp_matrix_raw)
        print()
        print(f"Numpy format {p}x{p} matrix:")
        print(syndomes_pxp_matrix)
        print()
        syndomes_pxp_matrix_sympy = Matrix(syndomes_pxp_matrix)
        print(f"Syndromes {p}x{p} sympy:")
        print(syndomes_pxp_matrix_sympy)
        print()
        try:
            inverse_syndromes_pxp_matrix_sympy = syndomes_pxp_matrix_sympy.inv_mod(251)
            inverse_syndromes_pxp_matrix = isaac_newton.matrix(inverse_syndromes_pxp_matrix_sympy)
            break
        except:
            print(">>> Failed... trying next p!")
            print()
            continue

    if inverse_syndromes_pxp_matrix is None:
        print("No p worked to create an invertible pxp matrix. ")
        exit(1)
    else:
        print(">>> Success!")
        print("Invertible matrix found for p = ", p, "; the inverted matrix is. ") # correct p now stored in this variable.
        print(inverse_syndromes_pxp_matrix)
        print()

    # finding ELP-coefficient matrix
    syndomes_px1_matrix_raw = []
    for i in range(p, 2*p):
        syndomes_px1_matrix_raw.append([-syndromes[i] % 251])
    syndromes_px1_matrix = isaac_newton.matrix(syndomes_px1_matrix_raw)
    print()
    print(f"Syndromes {p}x1 vector:")
    print(syndromes_px1_matrix)
    print()

    elp_coefficients_vector = isaac_newton.mod(inverse_syndromes_pxp_matrix * syndromes_px1_matrix, 251)
    print("ELP coefficients vector:")
    print(elp_coefficients_vector)
    print()

    elp_coefficients = list(elp_coefficients_vector)
    for i in range(len(elp_coefficients)):
        elp_coefficients[i] = int(elp_coefficients[i][0, 0])
    elp_coefficients.reverse()
    elp_coefficients.insert(0, 1)
    print("ELP coeffcients:", elp_coefficients)
    print()

    error_correcting_polynomial = polynomial.Polynomial(elp_coefficients)
    print("Error-correcting polynomial:")
    print(error_correcting_polynomial)
    print()
    return error_correcting_polynomial

def solve_elp(error_correcting_polynomial):
    Xr_inverses = []
    for i in range(251):
        # Get the raw float value
        raw_val = error_correcting_polynomial(i)

        # Check if the closest whole number modulo 251 equals 0
        if abs(round(raw_val) - raw_val) < 1e-7:
            if round(raw_val) % 251 == 0:
                Xr_inverses.append(i)

    print("Roots of the ELP, i.e., inverses of Xr's (error locations):", Xr_inverses)
    print(error_correcting_polynomial)
    print()
    fermat_inverse_vec = isaac_newton.vectorize(lambda a, p: pow(int(a), int(p - 2), int(p)))
    Xr = list(fermat_inverse_vec(Xr_inverses, 251))
    Xr = [int(i) for i in Xr]
    print("Hence, Xr are:", Xr)
    print()
    return Xr

def find_error_polynomial_coefficients_Yr(syndromes, Xr, p):
    augmented_matrix = []
    for row_power in range(1, p+1):
        row = []
        for i in Xr:
            i = i**row_power
            i = i % 251
            row.append(i)
        row.append(syndromes[row_power-1])
        augmented_matrix.append(row)
    print('Augmented matrix with Xr and syndromes vector')
    print(isaac_newton.matrix(augmented_matrix))
    print()
    # we've generated the raw LHS matrix with Xr's. Time to generate the matrix for syndromes on the RHS

    # great! time to do the scary bit. Gaussian elimination...
    # inverses: inverse_n = pow(n, -1, 251)

    for row_being_reduced in range(1, p + 1):
        pivot = augmented_matrix[row_being_reduced - 1][row_being_reduced - 1]
        print(pivot)  # when row=1, this should return the 0,0 element (top left)

        pivot_inverse = pow(pivot, -1, 251)
        print(pivot_inverse)
        for i in range(p + 1):  # multiplying the row by the inverse
            augmented_matrix[row_being_reduced - 1][i] *= pivot_inverse
            augmented_matrix[row_being_reduced - 1][i] %= 251
        print(f">>> Reduced row {row_being_reduced}. New matrix:")
        print(augmented_matrix)
        # now need to add to all other rows
        for row_to_add_to in range(1, p + 1):
            if row_being_reduced == row_to_add_to:
                continue
            else:
                print(">>> Mirror of pivot in target row:")
                # print(augmented_matrix[row_to_add_to - 1][row_being_reduced - 1])
                how_many_times_to_add = 251 - augmented_matrix[row_to_add_to - 1][row_being_reduced - 1]
                print(f">>> Adding row {row_being_reduced} to row {row_to_add_to} exactly {how_many_times_to_add} times.")
                for j in range(p + 1):
                    augmented_matrix[row_to_add_to - 1][j] += how_many_times_to_add * \
                                                              augmented_matrix[row_being_reduced - 1][j]
                    augmented_matrix[row_to_add_to - 1][j] %= 251
                print("Intermediate result:")
                print(isaac_newton.matrix(augmented_matrix))
                print()
    print(">>> Finished Gaussian Elimination! ")
    print(isaac_newton.matrix(augmented_matrix))
    print()

    Yr = []
    for i in range(len(Xr)):
        Yr.append(augmented_matrix[i][p])
    print(">>> Yr's found:")
    print(Yr)
    print()
    return Yr

def generate_error_polynomial_coefficients(Xr, Yr, n):
    ir = []
    for X in Xr:
        for i in range(250):
            test = pow(alpha, i, 251)
            if test == X:
                ir.append(i)

    print("Error locations (ir): ")
    print(ir)
    print()

    # Generating base list for coefficients
    error_polynomial_coefficients = []
    for i in range(n):
        error_polynomial_coefficients.append(0)

    for i in range(len(Xr)):
        error_polynomial_coefficients[ir[i]] = Yr[i]

    print("Error polynomial coefficients")
    print(error_polynomial_coefficients)
    print()

    return error_polynomial_coefficients


def decode(n, k, received):
    print()
    print("Len received:", len(received))
    print()

    if len(received) != n:
        print("Wrong received word length.")

    r_polynomial = polynomial.Polynomial(received)
    print("R_poly:", r_polynomial)
    print()

    syndromes = find_syndromes(n, k, r_polynomial)
    t = int(isaac_newton.floor((n-k)/2))
    error_correcting_polynomial = find_elp(syndromes, t)
    Xr = solve_elp(error_correcting_polynomial)
    Yr = find_error_polynomial_coefficients_Yr(syndromes, Xr, len(Xr))
    error_polynomial_coeffs = generate_error_polynomial_coefficients(Xr, Yr, n)
    print()
    print("Received polynomial coeffs:")
    received_polynomial_coeffs = list(int(i) for i in r_polynomial.coef[0])
    print(received_polynomial_coeffs)
    print()

    message_polynomial_coeffs = []
    for i in range(n):
        message_polynomial_coeffs.append((received_polynomial_coeffs[i] - error_polynomial_coeffs[i]) % 251)

    return message_polynomial_coeffs

from itertools import product

def reduce_poly_mod_251(poly):
    poly = list(int(i) for i in poly.coef)

    for i in range(len(poly)):
        poly[i] = poly[i] % 251

    # convert back into polynomial
    return polynomial.Polynomial(poly)

def find_message(decoded_coeff, n, k):
    generator_polynomial = polynomial.Polynomial([-int(alpha) % 251, 1])

    for i in range(2, n - k + 1):
        next_root_poly = polynomial.Polynomial([-pow(alpha, i, 251), 1])
        generator_polynomial = polymul_mod(generator_polynomial, next_root_poly, 251)

    print("Codeword polynomial")
    print(polynomial.Polynomial(decoded_coeff))
    print()
    print("Generator polynomial:")
    print(generator_polynomial)
    print()
    # Now need to divide decoded_coeff polynomial by generator_polynomial in F251...

    message_coeffs = []
    for i in range(n-k):
        degree_being_reduced = n-i
        coefficient_being_reduced = decoded_coeff[n-i-1]
        # since the leading coefficient of the generator polynomial is 1, the subsequent
        # coefficient of the quotient polynomial will simply be equal to coefficient_being_reduced.
        message_coeffs.insert(0, coefficient_being_reduced)

        # now need to actually subtract all the other crap.
        # first, find the quotient monomial, i.e., the current element of the quotient polynomial.
        # it will have degree n-1-i. Trust
        # thus, the list must go to n-2-i, as 0 is included, so the n-2-ith element will be the coeff.

        quotient_mono_coeffs = []
        for j in range(n-k-i):
            quotient_mono_coeffs.append(0)

        quotient_mono_coeffs[n-k-1-i] = coefficient_being_reduced
        quotient_mono = polynomial.Polynomial(quotient_mono_coeffs)
        print("Current quotient monomial: ")
        print(quotient_mono.trim())
        print()

        subtraction_polynomial = reduce_poly_mod_251(quotient_mono * generator_polynomial)

        decoded_poly = reduce_poly_mod_251(polynomial.Polynomial(decoded_coeff) - subtraction_polynomial)

        print("New (reduced) polynomial:")
        print(decoded_poly)
        decoded_coeff = list(int(i) for i in decoded_poly.coef)
        print(decoded_coeff)
        print()
        print("Message:")
        print(message_coeffs)





def run(n, k, message, e):
    codeword = encode(n, k, message)
    print()
    print("Codeword coefficients:")
    codeword_coeffs = list(int(i) for i in codeword.coef)
    print(codeword_coeffs)
    print()

    print(f"Received word, containing errors. ")
    transmitted = transmit(codeword, e)
    print(transmitted)
    print()

    decoded = decode(n, k, transmitted[0])
    print()
    print("Decoded polynomial coeffs:")
    print(decoded)
    print()
    print("Codeword coefficients reminder:")
    codeword_coeffs = list(int(i) for i in codeword.coef)
    print(codeword_coeffs)
    print()
    print("Error locations (real):")
    print(transmitted[1])
    print("Error values (real):")
    print(transmitted[2])

    if decoded == codeword_coeffs:
        print("Success! Codeword was error-corrected properly. ")
        print()
    else:
        print("I'm sorry bro we've messed this one up. ")
        exit(1)

    print("Moving onto decoding the original message... ")

    find_message(decoded, n, k)


def test():
    find_error_polynomial_coefficients_Yr([3, 6, 9], [2, 3], 2)

run(20, 10, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 3)
#test()



"""
augmented_matrix = [[1, 2, 3],
                    [4, 5, 6],
                    [7, 8, 9]]
def test_gaussian_elim_mod_251(n, k, p, augmented_matrix):
    for row_being_reduced in range(1, p + 1):
        pivot = augmented_matrix[row_being_reduced - 1][row_being_reduced - 1]
        print(pivot)# when row=1, this should return the 0,0 element (top left)

        pivot_inverse = pow(pivot, -1, 251)
        print(pivot_inverse)
        for i in range(p + 1):  # multiplying the row by the inverse
            augmented_matrix[row_being_reduced-1][i] *= pivot_inverse
            augmented_matrix[row_being_reduced-1][i] %= 251
        print(f"reduced row {row_being_reduced}. New matrix:")
        print(augmented_matrix)
        # now need to add to all other rows
        for row_to_add_to in range(1, p + 2):
            if row_being_reduced == row_to_add_to:
                continue
            else:
                print("mirror of pivot in target row:")
                print(augmented_matrix[row_to_add_to-1][row_being_reduced-1])
                how_many_times_to_add = 251 - augmented_matrix[row_to_add_to-1][row_being_reduced-1]
                print(f"adding row {row_being_reduced} to row {row_to_add_to} exactly {how_many_times_to_add} times.")
                for j in range(p + 1):
                    augmented_matrix[row_to_add_to-1][j] += how_many_times_to_add*augmented_matrix[row_being_reduced-1][j]
                    augmented_matrix[row_to_add_to-1][j] %= 251
                print(augmented_matrix)
                print()

    Yr = []
    for i in range(n - k):
        Yr.append(augmented_matrix[i][p])

    return Yr
print(test_gaussian_elim_mod_251(4, 1, 2, augmented_matrix))
"""
