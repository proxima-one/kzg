use core::borrow::Borrow;
use core::cmp::{Eq, PartialEq};
use core::ops::{Add, AddAssign, Div, Mul, Range};
use pairing::{
    group::{ff::Field, Curve, Group},
    Engine,
};

#[derive(Clone, Debug)]
pub struct Polynomial<E: Engine, const MAX_DEGREE: usize> {
    pub degree: usize,
    pub coeffs: [E::Fr; MAX_DEGREE],
}

impl<E: Engine, const MAX_DEGREE: usize> PartialEq<Polynomial<E, MAX_DEGREE>>
    for Polynomial<E, MAX_DEGREE>
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

impl<E: Engine, const MAX_DEGREE: usize> Eq for Polynomial<E, MAX_DEGREE> {}

pub struct PolynomialSlice<'a, E: Engine, const MAX_DEGREE: usize> {
    degree: usize,
    coeffs: &'a [E::Fr],
}

impl<E: Engine, const MAX_DEGREE: usize> Polynomial<E, MAX_DEGREE> {
    pub fn is_zero(&self) -> bool {
        self.degree() == 0 && self.coeffs[0] == E::Fr::zero()
    }

    pub fn new_zero() -> Polynomial<E, MAX_DEGREE> {
        Polynomial {
            degree: 0,
            coeffs: [E::Fr::zero(); MAX_DEGREE],
        }
    }

    pub fn new(coeffs: [E::Fr; MAX_DEGREE]) -> Polynomial<E, MAX_DEGREE> {
        // figure out what the initial degree is
        let degree = Self::compute_degree(&coeffs, MAX_DEGREE);
        Polynomial { degree, coeffs }
    }

    /// note: use this carefully, as setting the degree incorrect can lead to the degree being inconsistent
    pub fn new_from_coeffs(
        coeffs: [E::Fr; MAX_DEGREE],
        degree: usize,
    ) -> Polynomial<E, MAX_DEGREE> {
        Polynomial { degree, coeffs }
    }

    pub fn slice(&self, range: Range<usize>) -> PolynomialSlice<E, MAX_DEGREE> {
        PolynomialSlice {
            degree: Self::compute_degree(&self.coeffs, range.len()),
            coeffs: &self.coeffs[range],
        }
    }

    pub fn compute_degree(coeffs: &[E::Fr], upper_bound: usize) -> usize {
        let mut i = upper_bound;
        loop {
            if coeffs[i] != E::Fr::zero() {
                break i + 1;
            } else if i == 0 {
                break 0;
            };

            i -= 1;
        }
    }

    pub fn shrink_degree(&mut self) {
        let degree = Self::compute_degree(&self.coeffs, self.degree);
        self.degree = degree;
    }

    pub fn fixup_degree(&mut self) {
        let degree = Self::compute_degree(&self.coeffs, MAX_DEGREE);
        self.degree = degree;
    }

    pub fn lead(&self) -> E::Fr {
        self.coeffs[self.degree - 1]
    }

    pub fn constant(&self) -> E::Fr {
        self.coeffs[0]
    }

    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn eval(&self, x: E::Fr) -> E::Fr {
        let mut res = self.coeffs[0];
        let mut term = x;

        for &coeff in self.coeffs.iter().skip(1) {
            res += coeff * term;
            term *= term;
        }

        res
    }

    pub fn long_division(
        &self,
        divisor: &Self,
    ) -> (Polynomial<E, MAX_DEGREE>, Option<Polynomial<E, MAX_DEGREE>>) {
        if self.is_zero() {
            (Self::new_zero(), Some(self.clone()))
        } else if divisor.is_zero() {
            panic!("divisor must not be zero!")
        } else {
            let mut remainder = self.clone();
            let mut quotient = Polynomial::new_from_coeffs(
                [E::Fr::zero(); MAX_DEGREE],
                self.degree() - divisor.degree(),
            );

            // inverse guaranteed to succeed because divisor isn't 0.
            let lead_inverse = divisor.lead().invert().unwrap();
            while !remainder.is_zero() && remainder.degree() >= divisor.degree() {
                let factor = remainder.lead() * lead_inverse;
                let i = remainder.degree() - divisor.degree();
                quotient.coeffs[i] = factor;

                for (j, &coeff) in divisor.coeffs.iter().enumerate() {
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

impl<'a, E: Engine, const MAX_DEGREE: usize> Add for &'a Polynomial<E, MAX_DEGREE> {
    type Output = Polynomial<E, MAX_DEGREE>;

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

impl<E: Engine, R: Borrow<Polynomial<E, MAX_DEGREE>>, const MAX_DEGREE: usize> AddAssign<R>
    for Polynomial<E, MAX_DEGREE>
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
        let x: Polynomial<Bls12, 5> =
            Polynomial::new([3.into(), 0.into(), -Scalar::from(5), 0.into(), 3.into()]);
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
    }
}
