src: <https://almondine-song-c43.notion.site/Less-Intuitive-Homomorphisms-bdd10f0f47eb4a4c9607adc248b3e4bc>

# Less Intuitive Homomorphisms

Remember, a homomorphism from group A to group B exists **if and only if** a function Φ exists that takes an element from A and returns and element from B and Φ(x □ y) = Φ(x) ◇ Φ(y) for all x and y in A, where **□** is the binary operator of A and **◇** is the binary operator of B. **The existence of Φ is sufficient for the homomorphism to exist.** When solving these problems, keep in mind the function Φ might seem trivial or obvious in retrospect. Keep in mind that Φ needs to work for every possible x and y in A, but it does not need to map to every possible element in B.

Some of the problems ask for a homomorphism from group A to B and then the following problem asks for a homomorphism from B to A. This does not mean you should necessarily try to derive the inverse of Φ. The strategy you use to map from group A to group B in some cases will be completely different when you map from group B to group A. All of the above is true of monoids as well.

## Integers under addition to integer powers of 2 under multiplication

## The monoid of strings (including the empty string) under concatenation to non-negative integers under addition

## Positive real numbers greater than zero under multiplication to all real numbers under addition

- Hint
    
    Natural logarithm
    

## All real numbers under addition to positive real numbers greater than zero under multiplication

## The group of integers under addition to matrices of integers under addition

## The group of integer matrices under addition to integers under addition

## The monoid of 2x2 matrices under multiplication to the monoid of integers under multiplication

- Hint
    
    Matrix determinants
    

## The monoid of integers under multiplication to the monoid of 2x2 matrices under multiplication

- Hint
    
    Matrix determinants. There are several matrices for a given determinant, but you only need to pick one.
    

## The group of rational numbers (excluding rational numbers where the denominator is a multiple of `p`) to addition modulo `p`.

## The group of vectors under element-wise addition to polynomials under addition (hard)

- Hint
    
    Lagrange interpolation
    

## The group of integers under addition modulo n to an elliptic curve over a finite field with curve order n (hard)