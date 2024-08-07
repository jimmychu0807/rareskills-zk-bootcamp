// SPDX-License-Identifier: MIT
pragma solidity >= 0.8.23;

struct ECPoint {
    uint256 x;
    uint256 y;
}

struct G2Point {
    uint256 x1;
    uint256 x2;
    uint256 y1;
    uint256 y2;
}

contract Contract {
    // bn128 curve
    uint256 public constant curve_order =
        21888242871839275222246405745257275088548364400416034343698204186575808495617;
    uint256 public constant G1_x = 1;
    uint256 public constant G1_y = 2;
    uint256 public constant G2_x1 = 11559732032986387107991004021392285783925812861821192530917403151452391805634;
    uint256 public constant G2_x2 = 10857046999023057135944570762232829481370756359578518086990519993285655852781;
    uint256 public constant G2_y1 = 4082367875863433681332203403145435568316851327593401208105741076214120093531;
    uint256 public constant G2_y2 = 8495653923123431417604973247489272438418190587263600148770280649306958101930;

    function scalar_mul(ECPoint memory pt, uint256 s) public view returns (ECPoint memory) {
        uint256[3] memory input;
        input[0] = pt.x;
        input[1] = pt.y;
        input[2] = s;
        (bool ok, bytes memory res) = address(0x07).staticcall(abi.encode(input));
        require(ok);

        return abi.decode(res, (ECPoint));
    }

    function pt_add(ECPoint memory pt1, ECPoint memory pt2) public view returns (ECPoint memory) {
        uint256[4] memory input;
        input[0] = pt1.x;
        input[1] = pt1.y;
        input[2] = pt2.x;
        input[3] = pt2.y;
        (bool ok, bytes memory res) = address(0x06).staticcall(abi.encode(input));
        require(ok);

        return abi.decode(res, (ECPoint));
    }

    function pairing(
        ECPoint memory aG1,
        G2Point memory bG2,
        ECPoint memory cG1,
        G2Point memory dG2
    ) public view returns (bool) {
        uint256[12] memory input = [
            aG1.x, aG1.y,
            bG2.x2, bG2.x1, bG2.y2, bG2.y1,
            cG1.x, cG1.y,
            dG2.x2, dG2.x1, dG2.y2, dG2.y1
        ];
        (bool ok, bytes memory res) = address(0x08).staticcall(abi.encode(input));
        require(ok);

        return abi.decode(res, (bool));
    }

    function mod_exp(uint256 x, uint256 y, uint256 k) public view returns (uint256) {
        uint256[6] memory input = [32, 32, 32, x, y, k];
        (bool ok, bytes memory res) = address(0x05).staticcall(abi.encode(input));
        require(ok);

        return abi.decode(res, (uint256));
    }

    function rationalAdd1(
        ECPoint calldata A,
        ECPoint calldata B,
        uint256 num,
        uint256 den
    ) public view returns (bool verified) {
        // Set:
        //   aG1 = (A + B) * -1
        //   bG2 = G2
        //   cG1 = G1 * (num * pow(den, -1, curve_order))
        //   dG2 = G2
        // Then:
        //   aG1 * bG2 + cG1 * dG2 = 0
        //   pairing(aG1, bG2, cG1, dG2) -> true
        ECPoint memory G1 = ECPoint(G1_x, G1_y);
        G2Point memory G2 = G2Point(G2_x2, G2_x1, G2_y2, G2_y1);

        ECPoint memory g1_sum = pt_add(A, B);
        ECPoint memory aG1 = scalar_mul(g1_sum, curve_order - 1);
        G2Point memory bG2 = G2;

        uint256 val = mulmod(num, mod_exp(den, curve_order - 2, curve_order), curve_order);
        ECPoint memory cG1 = scalar_mul(G1, val);
        G2Point memory dG2 = G2;

        verified = pairing(aG1, bG2, cG1, dG2);
    }

    function rationalAdd2(
        ECPoint calldata A,
        ECPoint calldata B,
        uint256 num,
        uint256 den
    ) public view returns (bool verified) {
        ECPoint memory aG1 = pt_add(A, B);

        uint256 val = mulmod(num, mod_exp(den, curve_order - 2, curve_order), curve_order);
        ECPoint memory G1 = ECPoint(G1_x, G1_y);
        ECPoint memory bG1 = scalar_mul(G1, val);

        verified = aG1.x == bG1.x && aG1.y == bG1.y;
    }

    function matmul(uint256[] calldata matrix,
        uint256 n, // n x n for the matrix
        ECPoint[] calldata s, // n elements
        uint256[] calldata o // n elements
    ) public view returns (bool verified) {
        // revert if dimensions don't make sense or the matrices are empty
        require(n > 0, "n has to be larger than 0");
        require(matrix.length == n * n, "mismatched matrix size");
        require(s.length == n, "mismatched number of ECPoints");
        require(o.length == n, "mismatched number of uint256 points");

        ECPoint memory G1 = ECPoint(G1_x, G1_y);

        // return true if Ms == o elementwise. You need to do n equality checks.
        // If you're lazy, you can hardcode n to 3, but it is suggested that you do this with a for loop
        for (uint256 i = 0; i < n; i++) {
            // LHS calc
            ECPoint[] memory ec_pts = new ECPoint[](n);
            for (uint256 j = 0; j < n; j++) {
                ec_pts[j] = scalar_mul(s[j], matrix[i * n + j]);
            }
            ECPoint memory acc = ec_pts[0];
            for (uint j = 1; j < n; j++) {
                acc = pt_add(acc, ec_pts[j]);
            }

            // RHS calc
            ECPoint memory rhs = scalar_mul(G1, o[i]);
            verified = acc.x == rhs.x && acc.y == rhs.y;
            if (!verified) {
                // can immediately return if the above equality doesn't hold.
                return verified;
            }
        }
        verified = true;
    }
}
