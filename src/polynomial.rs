use core::borrow::Borrow;
use core::cmp::{Eq, PartialEq};
use core::ops::{Add, AddAssign, Div, Mul, Range};
use pairing::{
    group::{ff::Field, Curve, Group},
    Engine,
};
use std::iter::Iterator;

#[derive(Clone, Debug)]
pub struct Polynomial<E: Engine, const MAX_COEFFS: usize> {
    pub degree: usize,
    pub coeffs: [E::Fr; MAX_COEFFS],
}

impl<E: Engine, const MAX_COEFFS: usize> PartialEq<Polynomial<E, MAX_COEFFS>>
    for Polynomial<E, MAX_COEFFS>
{
    fn eq(&self, other: &Self) -> bool {
        if self.degree() != other.degree() {
            false
        } else {
            self.coeffs
                .iter()
                .zip(other.coeffs.iter())
                .all(|(l, r)| l == r)
        }
    }
}

impl<E: Engine, const MAX_COEFFS: usize> Eq for Polynomial<E, MAX_COEFFS> {}

pub struct PolynomialSlice<'a, E: Engine, const MAX_COEFFS: usize> {
    degree: usize,
    coeffs: &'a [E::Fr],
}

impl<E: Engine, const MAX_COEFFS: usize> Polynomial<E, MAX_COEFFS> {
    pub fn is_zero(&self) -> bool {
        self.degree() == 0 && self.coeffs[0] == E::Fr::zero()
    }

    pub fn new_zero() -> Polynomial<E, MAX_COEFFS> {
        Polynomial {
            degree: 0,
            coeffs: [E::Fr::zero(); MAX_COEFFS],
        }
    }

    pub fn new(coeffs: [E::Fr; MAX_COEFFS]) -> Polynomial<E, MAX_COEFFS> {
        // figure out what the initial degree is
        let degree = Self::compute_degree(&coeffs, MAX_COEFFS - 1);
        Polynomial { degree, coeffs }
    }

    /// note: use this carefully, as setting the degree incorrect can lead to the degree being inconsistent
    pub fn new_from_coeffs(
        coeffs: [E::Fr; MAX_COEFFS],
        degree: usize,
    ) -> Polynomial<E, MAX_COEFFS> {
        Polynomial { degree, coeffs }
    }

    pub fn slice(&self, range: Range<usize>) -> PolynomialSlice<E, MAX_COEFFS> {
        PolynomialSlice {
            degree: Self::compute_degree(&self.coeffs, range.len()),
            coeffs: &self.coeffs[range],
        }
    }

    pub fn compute_degree(coeffs: &[E::Fr], upper_bound: usize) -> usize {
        let mut i = upper_bound;
        loop {
            if i == 0 {
                break 0;
            } else if coeffs[i] != E::Fr::zero() {
                break i;
            }

            i -= 1;
        }
    }

    pub fn shrink_degree(&mut self) {
        let degree = Self::compute_degree(&self.coeffs, self.degree);
        self.degree = degree;
    }

    pub fn fixup_degree(&mut self) {
        let degree = Self::compute_degree(&self.coeffs, MAX_COEFFS - 1);
        self.degree = degree;
    }

    pub fn lead(&self) -> E::Fr {
        self.coeffs[self.degree]
    }

    pub fn constant(&self) -> E::Fr {
        self.coeffs[0]
    }

    pub fn num_coeffs(&self) -> usize {
        self.degree + 1
    }

    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn iter_coeffs(&self) -> impl Iterator<Item = &E::Fr> {
        self.coeffs.iter().take(self.num_coeffs())
    }

    pub fn eval(&self, x: E::Fr) -> E::Fr {
        let mut res = E::Fr::zero();
        let mut term = E::Fr::one();

        for &coeff in self.iter_coeffs() {
            res += coeff * term;
            term *= x;
        }

        res
    }

    pub fn long_division(
        &self,
        divisor: &Self,
    ) -> (
        Polynomial<E, MAX_COEFFS>,
        Option<Polynomial<E, MAX_COEFFS>>,
    ) {
        if self.is_zero() {
            (Self::new_zero(), Some(self.clone()))
        } else if divisor.is_zero() {
            panic!("divisor must not be zero!")
        } else {
            let mut remainder = self.clone();
            let mut quotient = Polynomial::new_from_coeffs(
                [E::Fr::zero(); MAX_COEFFS],
                self.degree() - divisor.degree(),
            );

            // inverse guaranteed to succeed because divisor isn't 0.
            let lead_inverse = divisor.lead().invert().unwrap();
            while !remainder.is_zero() && remainder.degree() >= divisor.degree() {
                let factor = remainder.lead() * lead_inverse;
                let i = remainder.degree() - divisor.degree();
                quotient.coeffs[i] = factor;

                for (j, &coeff) in divisor.iter_coeffs().enumerate() {
                    remainder.coeffs[i + j] -= coeff * factor;
                }

                remainder.shrink_degree();
            }

            quotient.fixup_degree();

            if remainder.is_zero() {
                (quotient, None)
            } else {
                (quotient, Some(remainder))
            }
        }
    }
}

impl<'a, E: Engine, const MAX_COEFFS: usize> Add for &'a Polynomial<E, MAX_COEFFS> {
    type Output = Polynomial<E, MAX_COEFFS>;

    fn add(self, rhs: Self) -> Self::Output {
        let (mut res, shorter) = if rhs.degree() > self.degree {
            (rhs.clone(), self)
        } else {
            (self.clone(), rhs)
        };

        for i in 0..shorter.degree() {
            res.coeffs[i] += shorter.coeffs[i];
        }

        res
    }
}

impl<E: Engine, R: Borrow<Polynomial<E, MAX_COEFFS>>, const MAX_COEFFS: usize> AddAssign<R>
    for Polynomial<E, MAX_COEFFS>
{
    fn add_assign(&mut self, rhs: R) {
        let rhs = rhs.borrow();
        for i in 0..rhs.degree() {
            self.coeffs[i] += rhs.coeffs[i];
        }

        if rhs.degree() > self.degree() {
            self.degree = rhs.degree();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bls12_381::{Bls12, Scalar};

    #[test]
    fn test_polynomial_division() {
        // test cases taken from https://tutorial.math.lamar.edu/Solutions/Alg/DividingPolynomials

        // 3x^4 - 5x^2 + 3 / x + 2 = 3x^3 - 6x^2 + 7x - 14 r 31
        let x: Polynomial<Bls12, 5> = Polynomial::new([
            3.into(),
            Scalar::zero(),
            -Scalar::from(5),
            Scalar::zero(),
            3.into(),
        ]);
        let y: Polynomial<Bls12, 5> = Polynomial::new([
            2.into(),
            Scalar::one(),
            Scalar::zero(),
            Scalar::zero(),
            Scalar::zero(),
        ]);

        let (q, r) = x.long_division(&y);
        assert!(r.is_some());
        assert_eq!(
            r.unwrap(),
            Polynomial::new([
                31.into(),
                Scalar::zero(),
                Scalar::zero(),
                Scalar::zero(),
                Scalar::zero()
            ])
        );
        assert_eq!(
            q,
            Polynomial::new([
                -Scalar::from(14),
                7.into(),
                -Scalar::from(6),
                3.into(),
                Scalar::zero(),
            ])
        );

        // x^3 + 2x^2 - 3x + 4 / x - 7 = x^2 + 9x + 60 r 424
        let x: Polynomial<Bls12, 4> =
            Polynomial::new([4.into(), -Scalar::from(3), 2.into(), Scalar::one()]);
        let y: Polynomial<Bls12, 4> = Polynomial::new([
            -Scalar::from(7),
            Scalar::one(),
            Scalar::zero(),
            Scalar::zero(),
        ]);

        let (q, r) = x.long_division(&y);
        assert!(r.is_some());
        assert_eq!(
            r.unwrap(),
            Polynomial::new([424.into(), Scalar::zero(), Scalar::zero(), Scalar::zero(),])
        );
        assert_eq!(
            q,
            Polynomial::new([60.into(), 9.into(), Scalar::one(), Scalar::zero(),])
        );

        // x^3 + 6x^2 + 13x + 10 / x + 2 = x^2 + 4x + 5 r 0
        let x: Polynomial<Bls12, 4> =
            Polynomial::new([10.into(), 13.into(), 6.into(), Scalar::one()]);
        let y: Polynomial<Bls12, 4> = Polynomial::new([
            Scalar::from(2),
            Scalar::one(),
            Scalar::zero(),
            Scalar::zero(),
        ]);

        let (q, r) = x.long_division(&y);
        assert!(r.is_none());
        assert_eq!(
            q,
            Polynomial::new([5.into(), 4.into(), Scalar::one(), Scalar::zero(),])
        );
    }

    #[test]
    fn test_eval_basic() {
        // y(x) = x^5 + 4x^3 + 7x^2 + 34
        let polynomial: Polynomial<Bls12, 6> = Polynomial::new([
            34.into(),
            Scalar::zero(),
            7.into(),
            4.into(),
            Scalar::zero(),
            Scalar::one(),
        ]);

        // y(0) = 34
        assert_eq!(polynomial.eval(Scalar::zero()), 34.into());
        // y(1) = 46
        assert_eq!(polynomial.eval(Scalar::one()), 46.into());
        // y(5) = 3834
        assert_eq!(polynomial.eval(5.into()), 3834.into());
    }
}
