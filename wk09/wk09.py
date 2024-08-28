import numpy as np
from functools import reduce
from py_ecc.bn128 import G1, G2, Z1, Z2, Z1, multiply, add, curve_order, pairing
import secrets
import galois

def trusted_setup(degree, GF, tau = None):
    if tau is None:
        # A rand between [1, ... GF_order - 1] inclusively
        tau = secrets.randbelow(GF.order - 1) + 1

    g1_pts = [multiply(G1, int((GF(tau) ** p))) for p in range(degree + 1)]
    g2_pts = [multiply(G2, int((GF(tau) ** p))) for p in range(degree + 1)]
    return (g1_pts, g2_pts)

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

    Returns:

    Raises:
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
        The return value is a polynomial with deg(poly) - 1
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

    return (Uw, Vw, Ww, H, T)

def main():
    # Define the Galois field
    print("Initializing galois field...")
    # For testing, switch to use a smaller field as below:
    # p = 59567
    p = curve_order
    GF = galois.GF(p)

    # Formula
    # out = 3x²y + 5xy - x - 2y + 3
    # witness = {x: 100, y: 100}
    
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
    x = to_galois(100, GF)
    y = to_galois(100, GF)
    v1 = to_galois(3, GF) * x * x
    v2 = v1 * y
    out = v2 + to_galois(5, GF) * x * y - x - to_galois(2, GF) * y + to_galois(3, GF)
    witness = GF(np.array([1, out, x, y, v1, v2]))

    # Convert ndarray into Galois field
    Lg = GF(np.array([[to_galois(x, GF) for x in row] for row in L], dtype="object"))
    Rg = GF(np.array([[to_galois(x, GF) for x in row] for row in R], dtype="object"))
    Og = GF(np.array([[to_galois(x, GF) for x in row] for row in O], dtype="object"))

    print("Converting r1cs to qap...")

    U, V, W, H, T = r1cs_to_qap(Lg, Rg, Og, witness, GF)
    HT = H * T

    print("Performing trusted setup...")

    # HT polynomial has the highest degree no. among (U, V, W, HT)
    g1_pp, g2_pp = trusted_setup(HT.degree, GF)

    print("Performing ecc computation on U, V, W, HT")

    Ug1 = ecc_eval(U, g1_pp, Z1)
    Vg2 = ecc_eval(V, g2_pp, Z2)
    Wg1 = ecc_eval(W, g1_pp, Z1)
    HTg1 = ecc_eval(HT, g1_pp, Z1)

    print("Performing 2 ecc pairings...")

    # Check if pairing (Vg2, Ug1) == pairing(G2, Wg1 + HTg1)
    lhs = pairing(Vg2, Ug1)
    rhs = pairing(G2, add(Wg1, HTg1))
    assert lhs == rhs, "pairing(Ug1, Vg2) != pairing(Wg1 + HTg1, G2)"

    print("pairing(Ug1, Vg2) == pairing(Wg1 + HTg1, G2)")

main()
