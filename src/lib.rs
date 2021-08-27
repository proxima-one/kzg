#[cfg(any(test, std))]
#[macro_use]
extern crate std;

use core::fmt::Debug;
use pairing::{
    group::{ff::Field, prime::PrimeCurveAffine, Curve},
    Engine,
};
use thiserror::Error;

pub mod polynomial;

use polynomial::{
    Polynomial,
    op_tree
};

/// parameters from tested setup
#[derive(Clone, Debug)]
pub struct KZGParams<E: Engine, const MAX_COEFFS: usize> {
    /// generator of g
    g: E::G1Affine,
    /// generator of G2
    h: E::G2Affine,
    /// g^alpha^1, g^alpha^2, ...
    gs: [E::G1Affine; MAX_COEFFS],
    /// g^alpha^1, g^alpha^2, ...
    hs: [E::G2Affine; MAX_COEFFS],
}

// the commitment - "C" in the paper. It's a single group element
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KZGCommitment<E: Engine>(E::G1Affine);
impl<E: Engine> Copy for KZGCommitment<E> {}

impl<E: Engine> KZGCommitment<E> {
    pub fn into_inner(self) -> E::G1Affine {
        self.0
    }

    pub fn inner(&self) -> &E::G1Affine {
        &self.0
    }
}

// A witness for a single element - "w_i" in the paper. It's a group element.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KZGWitness<E: Engine>(E::G1Affine);
impl<E: Engine> Copy for KZGWitness<E> {}

impl<E: Engine> KZGWitness<E> {
    pub fn into_inner(self) -> E::G1Affine {
        self.0
    }

    pub fn inner(&self) -> &E::G1Affine {
        &self.0
    }
}

// A witness for a several elements - "w_B" in the paper. It's a single group element plus a polynomial
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KZGBatchWitness<E: Engine, const MAX_COEFFS: usize>{
    r: Polynomial<E, MAX_COEFFS>,
    w: E::G1Affine
}

impl<E: Engine, const MAX_COEFFS: usize> KZGBatchWitness<E, MAX_COEFFS> {
    pub fn elem(self) -> E::G1Affine {
        self.w
    }

    pub fn elem_ref(&self) -> &E::G1Affine {
        &self.w
    }

    pub fn polynomial(&self) -> &Polynomial<E, MAX_COEFFS> {
        &self.r
    }
}

#[derive(Error, Debug)]
pub enum KZGError {
    #[error("no polynomial!")]
    NoPolynomial,
    #[error("point not on polynomial!")]
    PointNotOnPolynomial,
    #[error("batch opening remainder is zero!")]
    BatchOpeningZeroRemainder,
}

#[derive(Debug, Clone)]
pub struct KZGProver<'params, E: Engine, const MAX_COEFFS: usize> {
    parameters: &'params KZGParams<E, MAX_COEFFS>,
    polynomial: Option<Polynomial<E, MAX_COEFFS>>,
    commitment: Option<KZGCommitment<E>>,
}

#[derive(Debug, Clone)]
pub struct KZGVerifier<'params, E: Engine, const MAX_COEFFS: usize> {
    parameters: &'params KZGParams<E, MAX_COEFFS>,
}

impl<'params, E: Engine, const MAX_COEFFS: usize> KZGProver<'params, E, MAX_COEFFS> {
    /// initializes `polynomial` to zero polynomial
    pub fn new(parameters: &'params KZGParams<E, MAX_COEFFS>) -> Self {
        Self {
            parameters,
            polynomial: None,
            commitment: None,
        }
    }

    pub fn commit(&mut self, polynomial: Polynomial<E, MAX_COEFFS>) -> KZGCommitment<E> {
        let commitment = op_tree(
            polynomial.num_coeffs(),
            &|i| {
                if i == 0 {
                    self.parameters.g * polynomial.coeffs[i]
                } else {
                    self.parameters.gs[i - 1] * polynomial.coeffs[i]
                }
            },
            &|a, b| a + b
        ).to_affine();

        self.polynomial = Some(polynomial);
        let commitment = KZGCommitment(commitment);
        self.commitment = Some(commitment);
        commitment
    }

    pub fn open(&self) -> Result<Polynomial<E, MAX_COEFFS>, KZGError> {
        self.polynomial.clone().ok_or(KZGError::NoPolynomial)
    }

    pub fn commitment(&self) -> Option<KZGCommitment<E>> {
        self.commitment
    }

    pub fn commitment_ref(&self) -> Option<&KZGCommitment<E>> {
        self.commitment.as_ref()
    }

    pub fn has_commitment(&self) -> bool {
        self.commitment.is_some()
    }

    pub fn polynomial(&self) -> Option<&Polynomial<E, MAX_COEFFS>> {
        self.polynomial.as_ref()
    }

    pub fn has_polynomial(&self) -> bool {
        self.polynomial.is_some()
    }

    pub fn create_witness(&mut self, (x, y): (E::Fr, E::Fr)) -> Result<KZGWitness<E>, KZGError> {
        match self.polynomial {
            None => Err(KZGError::NoPolynomial),
            Some(ref polynomial) => {
                let mut dividend = polynomial.clone();
                dividend.coeffs[0] -= y;

                let mut divisor = Polynomial::new_from_coeffs([E::Fr::zero(); MAX_COEFFS], 1);
                divisor.coeffs[0] = -x;
                divisor.coeffs[1] = E::Fr::one();
                match dividend.long_division(&divisor) {
                    // by polynomial remainder theorem, if (x - point.x) does not divide self.polynomial, then
                    // self.polynomial(point.y) != point.1
                    (_, Some(_)) => Err(KZGError::PointNotOnPolynomial),
                    (psi, None) => {
                        let witness = op_tree(
                            psi.num_coeffs(), 
                            &|i| {
                                if i == 0 {
                                    self.parameters.g * psi.coeffs[i]
                                } else {
                                    self.parameters.gs[i - 1] * psi.coeffs[i]
                                }
                            },
                            &|a, b| a + b
                        ).to_affine(); 

                        Ok(KZGWitness(witness))
                    }
                }
            }
        }
    }

    // #[cfg(any(std))]
    pub fn create_witness_batched(&mut self, points: &Vec<(E::Fr, E::Fr)>) -> Result<KZGBatchWitness<E, MAX_COEFFS>, KZGError> {
        match self.polynomial {
            None => Err(KZGError::NoPolynomial),
            Some(ref polynomial) => {
                let zeros: Polynomial<E, MAX_COEFFS> = op_tree(
                    points.len(),
                    &|i| {
                        let mut coeffs = [E::Fr::zero(); MAX_COEFFS];
                        coeffs[0] = -points[i].0;
                        coeffs[1] = E::Fr::one();
                        Polynomial::new_from_coeffs(coeffs, 1)
                    },
                    &|a, b| a * b
                );
                
                let (xs, ys): (Vec<E::Fr>, Vec<E::Fr>) = points.iter().cloned().unzip();
                let interpolation = Polynomial::<E, MAX_COEFFS>::lagrange_interpolation(xs.as_slice(), ys.as_slice());

                let numerator = polynomial - &interpolation;
                let (psi, rem) = numerator.long_division(&zeros);
                match rem {
                    Some(_) => Err(KZGError::PointNotOnPolynomial),
                    None => {
                        let w = op_tree(
                            psi.num_coeffs(),
                            &|i| {
                                if i == 0 {
                                    self.parameters.g * psi.coeffs[i]
                                } else {
                                    self.parameters.gs[i - 1] * psi.coeffs[i]
                                }
                            },
                            &|a, b| a + b
                        ).to_affine();
                        Ok(KZGBatchWitness { r: interpolation, w })
                    }
                }
            }
        }
    }
}

impl<'params, E: Engine, const MAX_COEFFS: usize> KZGVerifier<'params, E, MAX_COEFFS> {
    pub fn new(parameters: &'params KZGParams<E, MAX_COEFFS>) -> Self {
        KZGVerifier { parameters }
    }

    pub fn verify_poly(
        &self,
        commitment: &KZGCommitment<E>,
        polynomial: &Polynomial<E, MAX_COEFFS>,
    ) -> bool {
        let check = op_tree(
            polynomial.num_coeffs(), 
            &|i| {
                if i == 0 {
                    self.parameters.g * polynomial.coeffs[i]
                } else {
                    self.parameters.gs[i - 1] * polynomial.coeffs[i]
                }
            },
            &|a, b| a + b
        );

        check.to_affine() == commitment.0
    }

    pub fn verify_eval(
        &self,
        (x, y): (E::Fr, E::Fr),
        commitment: &KZGCommitment<E>,
        witness: &KZGWitness<E>,
    ) -> bool {
        let lhs = E::pairing(
            &witness.0,
            &(self.parameters.hs[0].to_curve() - self.parameters.h * x).to_affine(),
        );
        let rhs = E::pairing(
            &(commitment.0.to_curve() - self.parameters.g * y).to_affine(),
            &self.parameters.h,
        );

        lhs == rhs
    }

    pub fn verify_eval_batched(
        &self,
        points: &Vec<(E::Fr, E::Fr)>,
        commitment: &KZGCommitment<E>,
        witness: &KZGBatchWitness<E, MAX_COEFFS>
    ) -> bool {
        let z: Polynomial<E, MAX_COEFFS> = op_tree(
            points.len(),
            &|i| {
                let mut coeffs = [E::Fr::zero(); MAX_COEFFS];
                coeffs[0] = -points[i].0;
                coeffs[1] = E::Fr::one();
                Polynomial::new_from_coeffs(coeffs, 1)
            },
            &|a, b| a * b
        );

        let hz = op_tree(
            z.num_coeffs(),
            &|i| {
                if i == 0 {
                    self.parameters.h * z.coeffs[i]
                } else {
                    self.parameters.hs[i - 1] * z.coeffs[i]
                }
            },
            &|a, b| a + b
        ).to_affine();

        let gr = op_tree(
            witness.r.num_coeffs(),
            &|i| {
                if i == 0 {
                    self.parameters.g * witness.r.coeffs[i]
                } else {
                    self.parameters.gs[i - 1] * witness.r.coeffs[i]
                }
            },
            &|a, b| a + b
        );

        let lhs = E::pairing(&witness.w, &hz);
        let rhs = E::pairing(
            &(commitment.0.to_curve() - gr).to_affine(),
            &self.parameters.h
        );

        lhs == rhs
    }
}

pub fn setup<E: Engine, const MAX_COEFFS: usize>(s: E::Fr) -> KZGParams<E, MAX_COEFFS> {
    let g = E::G1Affine::generator();
    let h = E::G2Affine::generator();

    let mut gs = [g; MAX_COEFFS];
    let mut hs = [h; MAX_COEFFS];

    let mut curr = g;

    for g in gs.iter_mut() {
        *g = (curr * s).to_affine();
        curr = *g;
    }

    let mut curr = h;
    for h in hs.iter_mut() {
        *h = (curr * s).to_affine();
        curr = *h;
    }

    KZGParams { g, h, gs, hs }
}

#[cfg(any(csprng_setup, test))]
use rand::random;

#[cfg(any(csprng_setup, test))]
pub fn csprng_setup<E: Engine, const MAX_COEFFS: usize>() -> KZGParams<E, MAX_COEFFS> {
    let s: E::Fr = random::<u64>().into();
    setup(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use bls12_381::{Bls12, Scalar};
    use lazy_static::lazy_static;
    use rand::{rngs::SmallRng, Rng, SeedableRng};
    use std::sync::Mutex;

    const RNG_SEED_0: [u8; 32] = [42; 32];
    const RNG_SEED_1: [u8; 32] = [79; 32];

    lazy_static! {
        static ref RNG_0: Mutex<SmallRng> = Mutex::new(SmallRng::from_seed(RNG_SEED_0));
        static ref RNG_1: Mutex<SmallRng> = Mutex::new(SmallRng::from_seed(RNG_SEED_1));
    }

    fn test_setup<E: Engine, const MAX_COEFFS: usize>() -> KZGParams<E, MAX_COEFFS> {
        let s: E::Fr = RNG_0.lock().unwrap().gen::<u64>().into();
        setup(s)
    }

    fn test_participants<'params, E: Engine, const MAX_COEFFS: usize>(params: &'params KZGParams<E, MAX_COEFFS>) -> (
        KZGProver<'params, E, MAX_COEFFS>, KZGVerifier<'params, E, MAX_COEFFS>
    ) {
        let prover = KZGProver::new(params);
        let verifier = KZGVerifier::new(params);

        (prover, verifier)
    }

    // never returns zero polynomial
    fn random_polynomial<E: Engine, const MAX_COEFFS: usize>(
        min_coeffs: usize,
    ) -> Polynomial<E, MAX_COEFFS> {
        let num_coeffs = RNG_1.lock().unwrap().gen_range(min_coeffs..MAX_COEFFS);
        let mut coeffs = [E::Fr::zero(); MAX_COEFFS];

        for i in 0..num_coeffs {
            coeffs[i] = RNG_1.lock().unwrap().gen::<u64>().into();
        }

        Polynomial::new_from_coeffs(coeffs, num_coeffs - 1)
    }

    fn assert_verify_poly<E: Engine + Debug, const MAX_COEFFS: usize>(
        verifier: &KZGVerifier<E, MAX_COEFFS>,
        commitment: &KZGCommitment<E>,
        polynomial: &Polynomial<E, MAX_COEFFS>,
    ) {
        assert!(
            verifier.verify_poly(&commitment, &polynomial),
            "verify_poly failed for commitment {:#?} and polynomial {:#?}",
            commitment,
            polynomial
        );
    }

    fn assert_verify_poly_fails<E: Engine + Debug, const MAX_COEFFS: usize>(
        verifier: &KZGVerifier<E, MAX_COEFFS>,
        commitment: &KZGCommitment<E>,
        polynomial: &Polynomial<E, MAX_COEFFS>,
    ) {
        assert!(
            !verifier.verify_poly(&commitment, &polynomial),
            "expected verify_poly to fail for commitment {:#?} and polynomial {:#?} but it didn't",
            commitment,
            polynomial
        );
    }

    fn assert_verify_eval<E: Engine + Debug, const MAX_COEFFS: usize>(
        verifier: &KZGVerifier<E, MAX_COEFFS>,
        point: (E::Fr, E::Fr),
        commitment: &KZGCommitment<E>,
        witness: &KZGWitness<E>,
    ) {
        assert!(
            verifier.verify_eval(point, &commitment, &witness),
            "verify_eval failed for point {:#?}, commitment {:#?}, and witness {:#?}",
            point,
            commitment,
            witness
        );
    }

    fn assert_verify_eval_fails<E: Engine + Debug, const MAX_COEFFS: usize>(
        verifier: &KZGVerifier<E, MAX_COEFFS>,
        point: (E::Fr, E::Fr),
        commitment: &KZGCommitment<E>,
        witness: &KZGWitness<E>,
    ) {
        assert!(!verifier.verify_eval(point, &commitment, &witness), "expected verify_eval to fail for for point {:#?}, commitment {:#?}, and witness {:#?}, but it didn't", point, commitment, witness);
    }

    #[test]
    fn test_basic() {
        let params = test_setup::<Bls12, 10>();
        let (mut prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(1);
        let commitment = prover.commit(polynomial.clone());

        assert_verify_poly(&verifier, &commitment, &polynomial);
        assert_verify_poly_fails(&verifier, &commitment, &random_polynomial(1));
    }

    fn random_field_elem_neq<E: Engine>(val: E::Fr) -> E::Fr {
        let mut v: E::Fr = RNG_1.lock().unwrap().gen::<u64>().into();
        while v == val {
            v = RNG_1.lock().unwrap().gen::<u64>().into();
        }

        v
    }

    #[test]
    fn test_modify_single_coeff() {
        let params = test_setup::<Bls12, 8>();
        let (mut prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(4);
        let commitment = prover.commit(polynomial.clone());

        let mut modified_polynomial = polynomial.clone();
        let new_coeff = random_field_elem_neq::<Bls12>(modified_polynomial.coeffs[2]);
        modified_polynomial.coeffs[2] = new_coeff;

        assert_verify_poly(&verifier, &commitment, &polynomial);
        assert_verify_poly_fails(&verifier, &commitment, &modified_polynomial);
    }

    #[test]
    fn test_eval_basic() {
        let params = test_setup::<Bls12, 13>();
        let (mut prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(5);
        let commitment = prover.commit(polynomial.clone());

        let x: Scalar = RNG_1.lock().unwrap().gen::<u64>().into();
        let y = polynomial.eval(x);

        let witness = prover.create_witness((x, y)).unwrap();
        assert_verify_eval(&verifier, (x, y), &commitment, &witness);

        let y_prime = random_field_elem_neq::<Bls12>(y);
        assert_verify_eval_fails(&verifier, (x, y_prime), &commitment, &witness);
    }

    #[test]
    fn test_eval_batched() {
        let params = test_setup::<Bls12, 15>();
        let (mut prover, verifier) = test_participants(&params);
        let polynomial = random_polynomial(8);
        let commitment = prover.commit(polynomial.clone());

        let mut points: Vec<(Scalar, Scalar)> = Vec::with_capacity(8);
        for _ in 0..8 {
            let x: Scalar = RNG_1.lock().unwrap().gen::<u64>().into();
            points.push((x, polynomial.eval(x)));
        }

        let witness = prover.create_witness_batched(&points).unwrap();
        assert!(verifier.verify_eval_batched(&points, &commitment, &witness));

        let mut other_points: Vec<(Scalar, Scalar)> = Vec::with_capacity(8);
        for _ in 0..8 {
            let x: Scalar = RNG_1.lock().unwrap().gen::<u64>().into();
            other_points.push((x, polynomial.eval(x)));
        }

        assert!(!verifier.verify_eval_batched(&other_points, &commitment, &witness))
    }
}
