from itertools import combinations
from sympy import symbols, interpolate, expand, simplify, Poly

evaluations = [5, 29, 47, 37, 29, 25, 29, 47]

def create_sets_of_6(evals):
    index_sets = list([list(c) for c in combinations(range(8), 6)])
    eval_sets = []
    for set in index_sets:
        eval_set = []
        for i in set:
            eval_set.append((i+1, evals[i]))
        eval_sets.append(eval_set)

    return eval_sets


points = create_sets_of_6(evaluations)
x = symbols('x')


for point_set in points:
    x = symbols('x')
    poly = expand(interpolate(point_set, x))
    degree = Poly(poly, x).degree()
    if degree == 3:
        P = Poly(poly, x)
        coeffs = P.all_coeffs()
        print(P)
        print(coeffs)

for i in [1, 2, 3, 4, 5, 6, 7, 8]:
    print(i**3 - 14*i**2 + 59*i - 41)