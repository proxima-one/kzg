use core::borrow::Borrow;
use core::cmp::{Eq, PartialEq};
use core::ops::{Add, AddAssign, Mul, Range, Sub, SubAssign};
use pairing::group::ff::{Field, PrimeField};
use std::collections::HashMap;
use std::iter::Iterator;

use crate::ft::EvaluationDomain;
use crate::utils::{log2_ceil, pad_to_power_of_two};
use crate::worker::Worker;

const FFT_MUL_THRESHOLD: usize = 128;

#[derive(Clone, Debug)]
pub struct Polynomial<S: PrimeField> {
    pub degree: usize,
    pub coeffs: Vec<S>,
}

impl<S: PrimeField> PartialEq<Polynomial<S>> for Polynomial<S> {
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

impl<S: PrimeField> Eq for Polynomial<S> {}

pub struct PolynomialSlice<'a, S: PrimeField> {
    pub degree: usize,
    pub coeffs: &'a [S],
}

impl<S: PrimeField> Polynomial<S> {
    pub fn is_zero(&self) -> bool {
        self.degree() == 0 && self.coeffs[0] == S::zero()
    }

    pub fn new_zero() -> Polynomial<S> {
        Polynomial {
            degree: 0,
            coeffs: vec![S::zero()],
        }
    }

    pub fn from_scalar(scalar: S) -> Polynomial<S> {
        Polynomial {
            degree: 0,
            coeffs: vec![scalar]
        }
    }

    pub fn new_monic_of_degree(degree: usize) -> Polynomial<S> {
        Polynomial {
            degree,
            coeffs: vec![S::one(); degree + 1]
        }
    }

    pub fn new_single_term(degree: usize) -> Polynomial<S> {
        let mut coeffs = vec![S::zero(); degree + 1];
        coeffs[degree] = S::one();
        Polynomial {
            degree,
            coeffs
        }
    }

    pub fn new_zero_with_size(cap: usize) -> Polynomial<S> {
        Polynomial {
            degree: 0,
            coeffs: vec![S::zero(); cap],
        }
    }

    pub fn new(coeffs: Vec<S>) -> Polynomial<S> {
        // figure out what the initial degree is
        let degree = Self::compute_degree(&coeffs, coeffs.len() - 1);
        Polynomial { degree, coeffs }
    }

    /// note: use this carefully, as setting the degree incorrect can lead to the degree being inconsistent
    pub fn new_from_coeffs(coeffs: Vec<S>, degree: usize) -> Polynomial<S> {
        Polynomial { degree, coeffs }
    }

    pub fn slice(&self, range: Range<usize>) -> PolynomialSlice<S> {
        PolynomialSlice {
            degree: Self::compute_degree(&self.coeffs, range.len()),
            coeffs: &self.coeffs[range],
        }
    }

    pub fn compute_degree(coeffs: &Vec<S>, upper_bound: usize) -> usize {
        let mut i = upper_bound;
        loop {
            if i == 0 {
                break 0;
            } else if coeffs[i] != S::zero() {
                break i;
            }

            i -= 1;
        }
    }

    pub fn truncate(&mut self, degree: usize) {
        self.degree = degree;
        self.coeffs.truncate(degree + 1);
    }

    pub fn reverse(&mut self) {
        self.coeffs.truncate(self.num_coeffs());
        self.coeffs.reverse();
    }

    pub fn shrink_degree(&mut self) {
        let degree = Self::compute_degree(&self.coeffs, self.degree);
        self.degree = degree;
    }

    pub fn fixup_degree(&mut self) {
        let degree = Self::compute_degree(&self.coeffs, self.coeffs.len() - 1);
        self.degree = degree;
    }

    pub fn lead(&self) -> S {
        self.coeffs[self.degree]
    }

    pub fn constant(&self) -> S {
        self.coeffs[0]
    }

    pub fn num_coeffs(&self) -> usize {
        self.degree + 1
    }

    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn iter_coeffs(&self) -> impl Iterator<Item = &S> {
        self.coeffs.iter().take(self.num_coeffs())
    }

    pub fn eval(&self, x: S) -> S {
        let mut res = self.coeffs[self.degree()];

        for i in (0..self.degree()).rev() {
            res *= x;
            res += self.coeffs[i];
        }

        res
    }

    pub fn fft_mul(&self, other: &Polynomial<S>, worker: &Worker) -> Polynomial<S> {
        let n = self.num_coeffs();
        let k = other.num_coeffs();
        let mut lhs = self.coeffs.clone();
        let mut rhs = other.coeffs.clone();
        lhs.resize(n + k, S::zero());
        rhs.resize(n + k, S::zero());

        let mut lhs = EvaluationDomain::from_coeffs(lhs).unwrap();
        let mut rhs = EvaluationDomain::from_coeffs(rhs).unwrap();

        lhs.fft(worker);
        rhs.fft(worker);
        lhs.mul_assign(worker, &rhs);
        lhs.ifft(worker);
        lhs.into()
    }

    pub fn best_mul(&self, other: &Polynomial<S>) -> Polynomial<S> {
        if self.degree() < FFT_MUL_THRESHOLD || other.degree() < FFT_MUL_THRESHOLD {
            self.clone() * other.clone()
        } else {
            let worker = Worker::new();
            self.fft_mul(&other, &worker)
        }
    }

    pub fn long_division(&self, divisor: &Self) -> (Polynomial<S>, Option<Polynomial<S>>) {
        if self.is_zero() {
            (Self::new_zero(), None)
        } else if divisor.is_zero() {
            panic!("divisor must not be zero!")
        } else if self.degree < divisor.degree() {
            (Self::new_zero(), Some(self.clone()))
        } else {
            let mut remainder = self.clone();
            let mut quotient = Polynomial::new_from_coeffs(
                vec![S::zero(); self.degree() - divisor.degree() + 1],
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

            if remainder.is_zero() {
                (quotient, None)
            } else {
                (quotient, Some(remainder))
            }
        }
    }

    pub fn fft_div(&self, divisor: &Self) -> (Polynomial<S>, Option<Polynomial<S>>) {
        let m = self.degree();
        let n = divisor.degree();

        let mut a_rev = self.clone();
        let mut b_rev = divisor.clone();

        a_rev.reverse();
        b_rev.reverse();

        let inv = b_rev.invert(m - n);
        
        let q_rev = a_rev.best_mul(&inv);
        let mut q = q_rev.clone();
        q.truncate(m - n);
        q.reverse();

        let r = self - &divisor.best_mul(&q);
        if r.is_zero() {
            (q, None)
        } else {
            (q, Some(r))
        }
    }

    /// computes the first degree + 1 terms of the formal series 1/f(x)
    /// panics if coeffs[0] == 0 or lead_coeff == 0
    pub fn invert(&self, degree: usize) -> Polynomial<S> {
        if degree == 0 {
            Polynomial::new_from_coeffs(vec![self.coeffs[0].invert().unwrap()], 0)
        } else {
            let c = self.invert(degree / 2);
            let mut res = c.best_mul(&Polynomial::from_scalar(2.into()).sub(&c.best_mul(self)));
            res.truncate(degree);
            res
        }
    }

    pub fn multi_eval(&self, xs: &[S]) -> Vec<S> {
        assert!(xs.len() > self.degree());
        let tree = SubProductTree::new_from_points(xs);
        tree.eval(xs.as_ref(), self)
    }

    pub fn lagrange_interpolation(xs: &[S], ys: &[S]) -> Polynomial<S> {
        assert_eq!(xs.len(), ys.len());

        if xs.len() == 1 {
            let coeffs = vec![ys[0] - xs[0], S::one()];
            return Polynomial::new_from_coeffs(coeffs, 1);
        }

        // let xs = pad_to_power_of_two(xs);
        // let ys = pad_to_power_of_two(ys);
        let tree = SubProductTree::new_from_points(xs);

        let mut m_prime = tree.product.clone();
        for i in 1..m_prime.num_coeffs() {
            m_prime.coeffs[i] *= S::from(i as u64);
        }
        m_prime.coeffs.remove(0);
        m_prime.degree -= 1;


        let cs: Vec<S> = m_prime.multi_eval(xs).iter().enumerate().map(|(i, c)| ys[i] * c.invert().unwrap()).collect();

        tree.linear_mod_combination(cs.as_slice())
    }

    pub fn scalar_multiplication(mut self, rhs: S) -> Polynomial<S> {
        for i in 0..self.num_coeffs() {
            self.coeffs[i] *= rhs;
        }
        self
    }
}

pub struct SubProductTree<S: PrimeField> {
    product: Polynomial<S>,
    left: Option<Box<SubProductTree<S>>>,
    right: Option<Box<SubProductTree<S>>>
}

impl<S: PrimeField> SubProductTree<S> {
    pub fn new_from_points(xs: &[S]) -> SubProductTree<S> {
        match xs.len() {
            1 => SubProductTree {
                product: Polynomial::new_from_coeffs(vec![-xs[0], S::one()], 1),
                left: None,
                right: None
            },
            n => {
                let left = SubProductTree::new_from_points(&xs[..n / 2]);
                let right = SubProductTree::new_from_points(&xs[n / 2..]);
                SubProductTree {
                    product: left.product.best_mul(&right.product),
                    left: Some(Box::new(left)),
                    right: Some(Box::new(right))
                }
            }
        }
    }

    pub fn eval(&self, xs: &[S], f: &Polynomial<S>) -> Vec<S> {
        let n = xs.len();

        if n == 1 {
            let y = f.eval(xs[0]);
            vec![y]
        } else {

            let left = self.left.as_ref().unwrap();
            let right = self.right.as_ref().unwrap();

            let (_, r0) = f.long_division(&left.product);
            let (_, r1) = f.long_division(&right.product);

            let mut l0 = left.eval(&xs[..n/2], &r0.unwrap());
            let l1 = right.eval(&xs[n/2..], &r1.unwrap());

            l0.extend(l1);
            l0
        }
    }

    pub fn linear_mod_combination(&self, cs: &[S]) -> Polynomial<S> {
        let n = cs.len();

        if n == 1 {
            Polynomial::new_from_coeffs(vec![cs[0]], 0)
        } else {
            let left = self.left.as_ref().unwrap();
            let right = self.right.as_ref().unwrap();

            let l = left.linear_mod_combination(&cs[..n/2]);
            let r = right.linear_mod_combination(&cs[n/2..]);
            
            right.product.best_mul(&l) + left.product.best_mul(&r)
        }
    }
}

fn op_tree_inner<T, F, O>(left: usize, size: usize, get_elem: &F, op: &O) -> T
where
    F: Fn(usize) -> T,
    O: Fn(T, T) -> T,
{
    assert!(size > 0);
    if size == 1 {
        get_elem(left)
    } else if size == 2 {
        op(get_elem(left), get_elem(left + 1))
    } else {
        let mid = left + (size / 2);
        op(
            op_tree_inner(left, size / 2, get_elem, op),
            op_tree_inner(mid, size - (size / 2), get_elem, op),
        )
    }
}

pub fn op_tree<T, F, O>(size: usize, get_elem: &F, op: &O) -> T
where
    F: Fn(usize) -> T,
    O: Fn(T, T) -> T,
{
    op_tree_inner(0, size, get_elem, op)
}

impl<'a, S: PrimeField> Add for &'a Polynomial<S> {
    type Output = Polynomial<S>;

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

impl<S: PrimeField> Add for Polynomial<S> {
    type Output = Polynomial<S>;

    fn add(self, rhs: Self) -> Self::Output {
        let (mut res, shorter) = if rhs.degree() > self.degree() {
            (rhs, self)
        } else {
            (self, rhs)
        };

        for i in 0..shorter.num_coeffs() {
            res.coeffs[i] += shorter.coeffs[i];
        }

        res
    }
}

impl<S: PrimeField, R: Borrow<Polynomial<S>>> AddAssign<R> for Polynomial<S> {
    fn add_assign(&mut self, rhs: R) {
        let rhs = rhs.borrow();
        for i in 0..rhs.num_coeffs() {
            self.coeffs[i] += rhs.coeffs[i];
        }

        if self.degree() < rhs.degree() {
            self.degree = rhs.degree();
        }
    }
}

impl<'a, S: PrimeField> Sub for &'a Polynomial<S> {
    type Output = Polynomial<S>;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut res = self.clone();
        if rhs.num_coeffs() > self.num_coeffs() {
            res.coeffs.resize(rhs.num_coeffs(), S::zero());
            res.degree = rhs.degree();
        }

        for i in 0..rhs.num_coeffs() {
            res.coeffs[i] -= rhs.coeffs[i];
        }

        res.shrink_degree();
        res
    }
}

impl<S: PrimeField, R: Borrow<Polynomial<S>>> SubAssign<R> for Polynomial<S> {
    fn sub_assign(&mut self, rhs: R) {
        let rhs = rhs.borrow();
        for i in 0..rhs.num_coeffs() {
            self.coeffs[i] -= rhs.coeffs[i];
        }

        self.fixup_degree()
    }
}

impl<S: PrimeField> Mul<Polynomial<S>> for Polynomial<S> {
    type Output = Polynomial<S>;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut res = Polynomial::new_zero_with_size(self.degree() + rhs.degree() + 1);
        for i in 0..self.num_coeffs() {
            for j in 0..rhs.num_coeffs() {
                res.coeffs[i + j] += self.coeffs[i] * rhs.coeffs[j];
            }
        }

        res.degree = self.degree() + rhs.degree();
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bls12_381::{Bls12, Scalar};

    #[test]
    fn test_long_division() {
        // test cases taken from https://tutorial.math.lamar.edu/Solutions/Alg/DividingPolynomials

        // 3x^4 - 5x^2 + 3 / x + 2 = 3x^3 - 6x^2 + 7x - 14 r 31
        let x = Polynomial::new(vec![
            3.into(),
            Scalar::zero(),
            -Scalar::from(5),
            Scalar::zero(),
            3.into(),
        ]);
        let y = Polynomial::new(vec![
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
            Polynomial::new(vec![
                31.into(),
                Scalar::zero(),
                Scalar::zero(),
                Scalar::zero(),
                Scalar::zero()
            ])
        );
        assert_eq!(
            q,
            Polynomial::new(vec![
                -Scalar::from(14),
                7.into(),
                -Scalar::from(6),
                3.into(),
                Scalar::zero(),
            ])
        );

        // x^3 + 2x^2 - 3x + 4 / x - 7 = x^2 + 9x + 60 r 424
        let x = Polynomial::new(vec![4.into(), -Scalar::from(3), 2.into(), Scalar::one()]);
        let y = Polynomial::new(vec![
            -Scalar::from(7),
            Scalar::one(),
            Scalar::zero(),
            Scalar::zero(),
        ]);

        let (q, r) = x.long_division(&y);
        assert!(r.is_some());
        assert_eq!(
            r.unwrap(),
            Polynomial::new(vec![
                424.into(),
                Scalar::zero(),
                Scalar::zero(),
                Scalar::zero(),
            ])
        );
        assert_eq!(
            q,
            Polynomial::new(vec![60.into(), 9.into(), Scalar::one(), Scalar::zero(),])
        );

        // x^3 + 6x^2 + 13x + 10 / x + 2 = x^2 + 4x + 5 r 0
        let x = Polynomial::new(vec![10.into(), 13.into(), 6.into(), Scalar::one()]);
        let y = Polynomial::new(vec![
            Scalar::from(2),
            Scalar::one(),
            Scalar::zero(),
            Scalar::zero(),
        ]);

        let (q, r) = x.long_division(&y);
        assert!(r.is_none());
        assert_eq!(
            q,
            Polynomial::new(vec![5.into(), 4.into(), Scalar::one(), Scalar::zero(),])
        );
    }

    #[test]
    fn test_fft_division() {
        // test cases taken from https://tutorial.math.lamar.edu/Solutions/Alg/DividingPolynomials

        // 3x^4 - 5x^2 + 3 / x + 2 = 3x^3 - 6x^2 + 7x - 14 r 31
        let x = Polynomial::new(vec![
            3.into(),
            Scalar::zero(),
            -Scalar::from(5),
            Scalar::zero(),
            3.into(),
        ]);
        let y = Polynomial::new(vec![
            2.into(),
            Scalar::one(),
            Scalar::zero(),
            Scalar::zero(),
            Scalar::zero(),
        ]);

        let (q, r) = x.fft_div(&y);
        assert!(r.is_some());
        assert_eq!(
            r.unwrap(),
            Polynomial::new(vec![
                31.into(),
                Scalar::zero(),
                Scalar::zero(),
                Scalar::zero(),
                Scalar::zero()
            ])
        );
        assert_eq!(
            q,
            Polynomial::new(vec![
                -Scalar::from(14),
                7.into(),
                -Scalar::from(6),
                3.into(),
                Scalar::zero(),
            ])
        );

        // x^3 + 2x^2 - 3x + 4 / x - 7 = x^2 + 9x + 60 r 424
        let x = Polynomial::new(vec![4.into(), -Scalar::from(3), 2.into(), Scalar::one()]);
        let y = Polynomial::new(vec![
            -Scalar::from(7),
            Scalar::one(),
            Scalar::zero(),
            Scalar::zero(),
        ]);

        let (q, r) = x.fft_div(&y);
        assert!(r.is_some());
        assert_eq!(
            r.unwrap(),
            Polynomial::new(vec![
                424.into(),
                Scalar::zero(),
                Scalar::zero(),
                Scalar::zero(),
            ])
        );
        assert_eq!(
            q,
            Polynomial::new(vec![60.into(), 9.into(), Scalar::one(), Scalar::zero(),])
        );

        // x^3 + 6x^2 + 13x + 10 / x + 2 = x^2 + 4x + 5 r 0
        let x = Polynomial::new(vec![10.into(), 13.into(), 6.into(), Scalar::one()]);
        let y = Polynomial::new(vec![
            Scalar::from(2),
            Scalar::one(),
            Scalar::zero(),
            Scalar::zero(),
        ]);

        let (q, r) = x.fft_div(&y);
        assert!(r.is_none());
        assert_eq!(
            q,
            Polynomial::new(vec![5.into(), 4.into(), Scalar::one(), Scalar::zero(),])
        );
    }

    #[test]
    fn test_eval_basic() {
        // y(x) = x^5 + 4x^3 + 7x^2 + 34
        let polynomial = Polynomial::new(vec![
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

    fn verify_tree(tree: &SubProductTree<Scalar>) {
        if tree.left.is_some() && tree.right.is_some() {
            assert!(
                tree.product == tree.left.as_ref().unwrap().product.best_mul(&tree.right.as_ref().unwrap().product)
            );
        }
    }

    #[test]
    fn test_new_subproduct_tree() {
        let xs = [Scalar::from(2), Scalar::from(5), Scalar::from(7), Scalar::from(90), Scalar::from(111), Scalar::from(31), Scalar::from(29)];

        let tree = SubProductTree::new_from_points(&xs);
        verify_tree(&tree);

        let xs = [Scalar::from(2), Scalar::from(5), Scalar::from(7), Scalar::from(90), Scalar::from(111)];
        let tree = SubProductTree::new_from_points(&xs);
        verify_tree(&tree);
    }

    #[test]
    fn test_fast_multi_eval() {
        let polynomial: Polynomial<Scalar> = Polynomial::new(
            vec![2, 5, 7, 90, 111]
                .into_iter()
                .map(|x| x.into())
                .collect(),
        );

        let xs: Vec<Scalar> = vec![1, 2, 3, 4, 5, 6, 7, 8].into_iter().map(|x| x.into()).collect();
        
        let mut fast = polynomial.multi_eval(xs.as_slice());
        fast.truncate(xs.len());
        let mut slow = Vec::new();

        for i in 0..xs.len() {
           let slow_y = polynomial.eval(xs[i]);
           slow.push(slow_y);
        }

        let slow: Vec<Scalar> =  xs.iter().map(|x| polynomial.eval(*x)).collect();
        assert!(fast == slow);
    }

    #[test]
    fn test_interpolation() {
        let xs: Vec<Scalar> = vec![2].into_iter().map(|x| x.into()).collect();
        let ys: Vec<Scalar> = vec![8].into_iter().map(|x| x.into()).collect();

        let interpolation = Polynomial::lagrange_interpolation(xs.as_slice(), ys.as_slice());

        for (&x, &y) in xs.iter().zip(ys.iter()) {
            assert_eq!(interpolation.eval(x), y);
        }

        let xs: Vec<Scalar> = vec![2, 5, 7, 90, 111, 31, 29]
            .into_iter()
            .map(|x| x.into())
            .collect();
        let ys: Vec<Scalar> = vec![8, 1, 43, 2, 87, 122, 13]
            .into_iter()
            .map(|x| x.into())
            .collect();
        let interpolation = Polynomial::lagrange_interpolation(xs.as_slice(), ys.as_slice());

        for (&x, &y) in xs.iter().zip(ys.iter()) {
            assert_eq!(interpolation.eval(x), y);
        }
    }
}
