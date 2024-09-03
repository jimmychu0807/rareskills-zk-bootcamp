import numpy as np
from functools import reduce
from py_ecc.bn128 import G1, G2, Z1, Z2, FQ12, multiply, add, neg, is_inf, curve_order, pairing, final_exponentiate
import secrets
import galois

def trusted_setup(degree, GF, U, V, W, tau = None, alpha = None, beta = None):
    """
    Generate the trusted setup parameters
    """
    if tau is None: tau = secrets.randbelow(GF.order - 1) + 1
    if alpha is None: alpha = secrets.randbelow(GF.order - 1) + 1
    if beta is None: beta = secrets.randbelow(GF.order - 1) + 1

    alpha_g1 = multiply(G1, alpha)
    beta_g2 = multiply(G2, beta)
    srs_g1 = [multiply(G1, int((GF(tau) ** p))) for p in range(degree + 1)]
    srs_g2 = [multiply(G2, int((GF(tau) ** p))) for p in range(degree + 1)]

    psi = [alpha * V[col] + beta * U[col] + W[col] for col in range(len(U))]
    psi_solved_g1 = [multiply(G1, int(p(tau))) for p in psi]

    return (alpha_g1, beta_g2, srs_g1, srs_g2, psi_solved_g1)

def ecc_eval(poly: galois.Poly, g_pp, zero_g):
    coeffs = poly.coeffs
    degrees = poly.degrees
    result = zero_g
    for idx, deg in enumerate(degrees):
        result = add(result, multiply(g_pp[deg], int(coeffs[idx])))
    return result

def to_galois(n: int, GF: galois.FieldArray):
    if n >= 0 and n < GF.order:
        return GF(n)
    elif n < 0:
        positive = int(n)
        while positive < 0: positive += GF.order
        return GF(positive)

    # n > GF.order
    positive = int(n)
    while positive >= GF.order: positive -= GF.order
    return GF(positive)

def print_polys(poly_arr, prefix = None):
    for idx, poly in enumerate(poly_arr):
        print(f"{prefix}[{idx}]: {poly}")


def r1cs_to_qap(L: np.array, R: np.array, O: np.array, w: np.array, GF: galois.FieldArray) -> (np.array, np.array, np.array, np.array, np.array):
    """
    Convert an R1CS to QAP. All the computation is done in finite field
    They need to satisfy this relationship: L·w ⊙ R·w = O·w

    Paramters:
    L (np.array): m x n array, wrapped in Galois field
    R (np.array): m x n array, wrapped in Galois field
    O (np.array): m x n array, wrapped in Galois field
    w (np.array): n x 1 vector, wrapped in Galois field
    GF (galois.FieldArray): The finite field order
    """

    # Check for the params size
    if (L.ndim != 2 or R.ndim != 2 or O.ndim != 2 or w.ndim != 1
        or L.shape[0] != R.shape[0] or R.shape[0] != O.shape[0]
        or L.shape[1] != R.shape[1] or R.shape[1] != O.shape[1]
        or w.shape[0] != O.shape[1]
       ):
        raise Exception("Parameters sizes mismatched.")

    def interpolate_col(col):
        xs = GF(np.array(range(1, len(col) + 1)))
        return galois.lagrange_poly(xs, col)

    def inner_product_polys_with_w(polys, witness):
        mul_ = lambda x, y: x * y
        sum_ = lambda x, y: x + y
        return reduce(sum_, map(mul_, polys, witness))

    def get_t_poly(poly: galois.Poly) -> galois.Poly:
        """
        Returns:
        The return value is a polynomial with deg(poly) - 1. We want to leave one degree of freedom
        for h(x).
        """
        return galois.Poly.Roots(list(range(1, poly.degree)), field = GF)

    Lg = GF(L)
    Rg = GF(R)
    Og = GF(O)

    U = np.apply_along_axis(interpolate_col, 0, Lg)
    V = np.apply_along_axis(interpolate_col, 0, Rg)
    W = np.apply_along_axis(interpolate_col, 0, Og)

    # Uw·Vw = Ww + h(x)t(x)
    Uw = inner_product_polys_with_w(U, w)
    Vw = inner_product_polys_with_w(V, w)
    Ww = inner_product_polys_with_w(W, w)

    # Uw·Vw - Ww
    result_poly = Uw * Vw - Ww
    T = get_t_poly(result_poly)
    H = result_poly // T

    return (Uw, Vw, Ww, U, V, W, H, T)

def sum_of_product(g1s, w):
    if (len(g1s) != len(w)): raise Exception("sum_of_product parameter sizes not equal.")

    mul_ = lambda i: multiply(g1s[i], int(w[i]))
    sum_ = lambda x, y: add(x, y)
    result = reduce(sum_, map(mul_, [i for i in range(len(g1s))]))

    return result

def main():
    # Define the Galois field
    print("Initializing galois field...")
    GF = galois.GF(curve_order)

    # Formula
    # out = 3x²y + 5xy - x - 2y + 3
    secret = { "x": 100, "y": 100 }

    # Define the matrices
    L = np.array([[0,0,3,0,0,0],
                  [0,0,0,0,1,0],
                  [0,0,1,0,0,0]], dtype="object")

    R = np.array([[0,0,1,0,0,0],
                  [0,0,0,1,0,0],
                  [0,0,0,5,0,0]], dtype="object")

    O = np.array([[0,0,0,0,1,0],
                  [0,0,0,0,0,1],
                  [-3,1,1,2,0,-1]], dtype="object")

    print("Computing witness and Lg, Rg, Og...")

    # Define the witness
    x = to_galois(secret["x"], GF)
    y = to_galois(secret["y"], GF)
    v1 = to_galois(3, GF) * x * x
    v2 = v1 * y
    out = v2 + to_galois(5, GF) * x * y - x - to_galois(2, GF) * y + to_galois(3, GF)
    witness = GF(np.array([1, out, x, y, v1, v2]))

    # Convert ndarray into Galois field
    Lg = GF(np.array([[to_galois(x, GF) for x in row] for row in L], dtype="object"))
    Rg = GF(np.array([[to_galois(x, GF) for x in row] for row in R], dtype="object"))
    Og = GF(np.array([[to_galois(x, GF) for x in row] for row in O], dtype="object"))

    print("Converting r1cs to qap...")

    Uw, Vw, Ww, U, V, W, H, T = r1cs_to_qap(Lg, Rg, Og, witness, GF)
    HT = H * T

    print("Performing trusted setup...")

    # HT polynomial has the highest degree no. among (U, V, W, HT)
    alpha_g1, beta_g2, srs_g1, srs_g2, psi_g1 = trusted_setup(HT.degree, GF, U, V, W)

    print("Performing prover steps (ecc computation)...")

    A_g1 = add(alpha_g1, ecc_eval(Uw, srs_g1, Z1))
    B_g2 = add(beta_g2, ecc_eval(Vw, srs_g2, Z2))

    HT_g1 = ecc_eval(HT, srs_g1, Z1)
    C_g1 = add(sum_of_product(psi_g1, witness), HT_g1)

    print("Performing verifier steps...")

    # Check if I₁₂ == neg([A]₁)·[B]₂ + [α]₁·[β]₂ + [C]₁·G₂
    A_B = pairing(B_g2, neg(A_g1))
    alpha_beta = pairing(beta_g2, alpha_g1)
    C_G2 = pairing(G2, C_g1)

    # Performing F12 field arithmetic: https://ethereum.stackexchange.com/questions/158662/how-to-add-fq12-points-result-of-pairings-when-using-ethereums-py-ecc-libra
    print("I₁₂ == neg([A]₁)·[B]₂ + [α]₁·[β]₂ + [C]₁·G₂:", FQ12.one() == final_exponentiate(A_B * alpha_beta * C_G2))

main()
