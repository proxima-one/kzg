use std::fmt::Debug;
use blstrs::{G1Affine, G2Affine, Scalar, pairing};
use pairing::group::{Curve, ff::Field, prime::PrimeCurveAffine};
use thiserror::Error;

#[cfg(feature = "serde_support")]
use serde::{Serialize, Deserialize};

pub mod utils;
pub mod ft;
pub mod polynomial;

use polynomial::{Polynomial, SubProductTree, op_tree};

/// parameters from tested setup
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde_support", derive(Serialize, Deserialize))]
pub struct KZGParams {
    /// generator of g
    g: G1Affine,
    /// generator of G2
    h: G2Affine,
    /// g^alpha^1, g^alpha^2, ...
    gs: Vec<G1Affine>,
    /// g^alpha^1, g^alpha^2, ...
    hs: Vec<G2Affine>,
}

/// the commitment - "C" in the paper. It's a single group element
pub type KZGCommitment = G1Affine;
/// A witness for a single element - "w_i" in the paper. It's a group element.
pub type KZGWitness = G1Affine;

// A witness for a several elements - "w_B" in the paper. It's a single group element plus a polynomial
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KZGBatchWitness {
    r: Polynomial,
    w: G1Affine,
}

impl KZGBatchWitness {
    pub fn elem(&self) -> G1Affine {
        self.w
    }

    pub fn elem_ref(&self) -> &G1Affine {
        &self.w
    }

    pub fn polynomial(&self) -> &Polynomial {
        &self.r
    }

    pub fn from_inner(r: Polynomial, w: G1Affine) -> Self {
        KZGBatchWitness { r, w }
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
    #[error("polynomial degree too large")]
    PolynomialDegreeTooLarge,
}

#[derive(Debug, Clone)]
pub struct KZGProver<'params> {
    parameters: &'params KZGParams,
    polynomial: Option<Polynomial>,
    commitment: Option<KZGCommitment>,
}

#[derive(Debug, Clone)]
pub struct KZGVerifier<'params> {
    parameters: &'params KZGParams,
}

impl<'params> KZGProver<'params> {
    /// initializes `polynomial` to zero polynomial
    pub fn new(parameters: &'params KZGParams) -> Self {
        Self {
            parameters,
            polynomial: None,
            commitment: None,
        }
    }

    pub fn parameters(&self) -> &'params KZGParams {
        self.parameters
    }

    pub fn commit(&mut self, polynomial: Polynomial) -> KZGCommitment {
        let mut commitment = self.parameters.g * polynomial.coeffs[0];
        for i in 0..polynomial.degree() {
            commitment += self.parameters.gs[i] * polynomial.coeffs[i + 1];
        }

        self.polynomial = Some(polynomial);
        let commitment = commitment.to_affine();
        self.commitment = Some(commitment);
        commitment
    }

    pub fn open(&self) -> Result<Polynomial, KZGError> {
        self.polynomial.clone().ok_or(KZGError::NoPolynomial)
    }

    pub fn commitment(&self) -> Option<KZGCommitment> {
        self.commitment
    }

    pub fn commitment_ref(&self) -> Option<&KZGCommitment> {
        self.commitment.as_ref()
    }

    pub fn set_commitment(&mut self, commitment: KZGCommitment, polynomial: Polynomial) {
        self.commitment = Some(commitment);
        self.polynomial = Some(polynomial);
    }

    pub fn has_commitment(&self) -> bool {
        self.commitment.is_some()
    }

    pub fn polynomial(&self) -> Option<&Polynomial> {
        self.polynomial.as_ref()
    }

    pub fn has_polynomial(&self) -> bool {
        self.polynomial.is_some()
    }

    pub fn create_witness(&self, (x, y): (Scalar, Scalar)) -> Result<KZGWitness, KZGError> {
        match self.polynomial {
            None => Err(KZGError::NoPolynomial),
            Some(ref polynomial) => {
                let mut dividend = polynomial.clone();
                dividend.coeffs[0] -= y;

                let divisor = Polynomial::new_from_coeffs(vec![-x, Scalar::one()], 1);
                match dividend.long_division(&divisor) {
                    // by polynomial remainder theorem, if (x - point.x) does not divide self.polynomial, then
                    // self.polynomial(point.y) != point.1
                    (_, Some(_)) => Err(KZGError::PointNotOnPolynomial),
                    (psi, None) => {
                        let mut w = self.parameters.g * psi.coeffs[0];
                        for i in 0..psi.degree() {
                            w += self.parameters.gs[i] * psi.coeffs[i + 1];
                        }

                        Ok(w.to_affine())
                    }
                }
            }
        }
    }

    pub fn create_witness_batched(
        &self,
        xs: &[Scalar],
        ys: &[Scalar]
    ) -> Result<KZGBatchWitness, KZGError> {
        match self.polynomial {
            None => Err(KZGError::NoPolynomial),
            Some(ref polynomial) => {
                let tree = SubProductTree::new_from_points(xs);

                let interpolation =
                    Polynomial::lagrange_interpolation_with_tree(xs, ys, &tree);
               
                let numerator = polynomial - &interpolation;
                let (psi, rem) = numerator.long_division(&tree.product);
                match rem {
                    Some(_) => Err(KZGError::PointNotOnPolynomial),
                    None => {
                        let mut w = self.parameters.g * psi.coeffs[0];
                        for i in 0..psi.degree() {
                            w += self.parameters.gs[i] * psi.coeffs[i + 1];
                        }
                        Ok(KZGBatchWitness {
                            r: interpolation,
                            w: w.to_affine(),
                        })
                    }
                }
            }
        }
    }
}

impl<'params> KZGVerifier<'params> {
    pub fn new(parameters: &'params KZGParams) -> Self {
        KZGVerifier { parameters }
    }

    pub fn verify_poly(
        &self,
        commitment: &KZGCommitment,
        polynomial: &Polynomial,
    ) -> bool {
        let mut check = self.parameters.g * polynomial.coeffs[0];
        for i in 0..polynomial.degree() {
            check += self.parameters.gs[i]* polynomial.coeffs[i + 1];
        }

        check.to_affine() == *commitment
    }

    pub fn verify_eval(
        &self,
        (x, y): (Scalar, Scalar),
        commitment: &KZGCommitment,
        witness: &KZGWitness,
    ) -> bool {
        let lhs = pairing(
            witness,
            &(self.parameters.hs[0].to_curve() - self.parameters.h * x).to_affine(),
        );
        let rhs = pairing(
            &(commitment.to_curve() - self.parameters.g * y).to_affine(),
            &self.parameters.h,
        );

        lhs == rhs
    }

    pub fn verify_eval_batched(
        &self,
        xs: &[Scalar],
        commitment: &KZGCommitment,
        witness: &KZGBatchWitness,
    ) -> bool {
        let z: Polynomial = op_tree(
            xs.len(),
            &|i| {
                let mut coeffs = vec![-xs[i], Scalar::one()];
                coeffs[0] = -xs[i];
                coeffs[1] = Scalar::one();
                Polynomial::new_from_coeffs(coeffs, 1)
            },
            &|a, b| a.best_mul(&b),
        );

        let mut hz = self.parameters.h * z.coeffs[0];
        for i in 0..z.degree() {
            hz += self.parameters.hs[i] * z.coeffs[i + 1];
        }

        let mut gr = self.parameters.g * witness.r.coeffs[0];
        for i in 0..witness.r.degree() {
            gr += self.parameters.gs[i] * witness.r.coeffs[i + 1];
        }

        let lhs = pairing(&witness.w, &hz.to_affine());
        let rhs = pairing(
            &(commitment.to_curve() - gr).to_affine(),
            &self.parameters.h,
        );

        lhs == rhs
    }
}

pub fn setup(s: Scalar, num_coeffs: usize) -> KZGParams {
    let g = G1Affine::generator();
    let h = G2Affine::generator();

    let mut gs = vec![g; num_coeffs];
    let mut hs = vec![h; num_coeffs];

    let mut curr = g;

    for g in gs.iter_mut() {
        *g =  (curr * s).to_affine();
        curr = *g;
    }

    let mut curr = h;
    for h in hs.iter_mut() {
        *h =  (curr * s).to_affine();
        curr = *h;
    }

    KZGParams { g, h, gs, hs }
}

#[cfg(any(csprng_setup, test))]
use rand::random;

#[cfg(any(csprng_setup, test))]
pub fn csprng_setup(num_coeffs: usize) -> KZGParams {
    let s: Scalar = random::<u64>().into();
    setup(s, num_coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use lazy_static::lazy_static;
    use rand::{rngs::SmallRng, Rng, SeedableRng};
    use std::sync::Mutex;

    const RNG_SEED_0: [u8; 32] = [42; 32];
    const RNG_SEED_1: [u8; 32] = [79; 32];

    lazy_static! {
        static ref RNG_0: Mutex<SmallRng> = Mutex::new(SmallRng::from_seed(RNG_SEED_0));
        static ref RNG_1: Mutex<SmallRng> = Mutex::new(SmallRng::from_seed(RNG_SEED_1));
    }

    fn test_setup<const MAX_COEFFS: usize>() -> KZGParams {
        let s: Scalar = RNG_0.lock().unwrap().gen::<u64>().into();
        setup(s, MAX_COEFFS)
    }

    fn test_participants<'params>(
        params: &'params KZGParams,
    ) -> (KZGProver<'params>, KZGVerifier<'params>) {
        let prover = KZGProver::new(params);
        let verifier = KZGVerifier::new(params);

        (prover, verifier)
    }

    // never returns zero polynomial
    fn random_polynomial(min_coeffs: usize, max_coeffs: usize) -> Polynomial {
        let num_coeffs = RNG_1.lock().unwrap().gen_range(min_coeffs..max_coeffs);
        let mut coeffs = vec![Scalar::zero(); max_coeffs];

        for i in 0..num_coeffs {
            coeffs[i] = RNG_1.lock().unwrap().gen::<u64>().into();
        }

        let mut poly = Polynomial::new_from_coeffs(coeffs, num_coeffs - 1);
        poly.shrink_degree();
        poly
    }

    fn assert_verify_poly(
        verifier: &KZGVerifier,
        commitment: &KZGCommitment,
        polynomial: &Polynomial,
    ) {
        assert!(
            verifier.verify_poly(&commitment, &polynomial),
            "verify_poly failed for commitment {:#?} and polynomial {:#?}",
            commitment,
            polynomial
        );
    }

    fn assert_verify_poly_fails(
        verifier: &KZGVerifier,
        commitment: &KZGCommitment,
        polynomial: &Polynomial,
    ) {
        assert!(
            !verifier.verify_poly(&commitment, &polynomial),
            "expected verify_poly to fail for commitment {:#?} and polynomial {:#?} but it didn't",
            commitment,
            polynomial
        );
    }

    fn assert_verify_eval(
        verifier: &KZGVerifier,
        point: (Scalar, Scalar),
        commitment: &KZGCommitment,
        witness: &KZGWitness,
    ) {
        assert!(
            verifier.verify_eval(point, &commitment, &witness),
            "verify_eval failed for point {:#?}, commitment {:#?}, and witness {:#?}",
            point,
            commitment,
            witness
        );
    }

    fn assert_verify_eval_fails(
        verifier: &KZGVerifier,
        point: (Scalar, Scalar),
        commitment: &KZGCommitment,
        witness: &KZGWitness,
    ) {
        assert!(!verifier.verify_eval(point, &commitment, &witness), "expected verify_eval to fail for for point {:#?}, commitment {:#?}, and witness {:#?}, but it didn't", point, commitment, witness);
    }

    #[test]
    fn test_basic() {
        let params = test_setup::<10>();
        let (mut prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(1, 10);
        let commitment = prover.commit(polynomial.clone());

        assert_verify_poly(&verifier, &commitment, &polynomial);
        assert_verify_poly_fails(&verifier, &commitment, &random_polynomial(1, 10));
    }

    fn random_field_elem_neq(val: Scalar) -> Scalar {
        let mut v: Scalar = RNG_1.lock().unwrap().gen::<u64>().into();
        while v == val {
            v = RNG_1.lock().unwrap().gen::<u64>().into();
        }

        v
    }

    #[test]
    fn test_modify_single_coeff() {
        let params = test_setup::<8>();
        let (mut prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(4, 8);
        let commitment = prover.commit(polynomial.clone());

        let mut modified_polynomial = polynomial.clone();
        let new_coeff = random_field_elem_neq(modified_polynomial.coeffs[2]);
        modified_polynomial.coeffs[2] = new_coeff;

        assert_verify_poly(&verifier, &commitment, &polynomial);
        assert_verify_poly_fails(&verifier, &commitment, &modified_polynomial);
    }

    #[test]
    fn test_eval_basic() {
        let params = test_setup::<13>();
        let (mut prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(5, 13);
        let commitment = prover.commit(polynomial.clone());

        let x: Scalar = RNG_1.lock().unwrap().gen::<u64>().into();
        let y = polynomial.eval(x);

        let witness = prover.create_witness((x, y)).unwrap();
        assert_verify_eval(&verifier, (x, y), &commitment, &witness);

        let y_prime = random_field_elem_neq(y);
        assert_verify_eval_fails(&verifier, (x, y_prime), &commitment, &witness);

        // test degree 1 edge case
        let mut coeffs = vec![Scalar::zero(); 13];
        coeffs[0] = 3.into();
        coeffs[1] = 1.into();
        let polynomial = Polynomial::new(coeffs);

        let commitment = prover.commit(polynomial);
        let witness = prover.create_witness((1.into(), 4.into())).unwrap();
        assert_verify_eval(&verifier, (1.into(), 4.into()), &commitment, &witness);
        assert_verify_eval_fails(&verifier, (1.into(), 5.into()), &commitment, &witness);
    }

    #[test]
    fn test_eval_batched() {
        let params = test_setup::<15>();
        let (mut prover, verifier) = test_participants(&params);
        let polynomial = random_polynomial(8, 15);
        let commitment = prover.commit(polynomial.clone());

        let mut xs: Vec<Scalar> = Vec::with_capacity(8);
        let mut ys: Vec<Scalar> = Vec::with_capacity(8);
        for _ in 0..8 {
            let x: Scalar = RNG_1.lock().unwrap().gen::<u64>().into();
            xs.push(x);
            ys.push(polynomial.eval(x));
        }

        let witness = prover.create_witness_batched(xs.as_slice(), ys.as_slice()).unwrap();
        assert!(verifier.verify_eval_batched(xs.as_slice(), &commitment, &witness));

        let mut xs: Vec<Scalar> = Vec::with_capacity(8);
        let mut ys: Vec<Scalar> = Vec::with_capacity(8);
        for _ in 0..8 {
            let x: Scalar = RNG_1.lock().unwrap().gen::<u64>().into();
            xs.push(x);
            ys.push(polynomial.eval(x));
        }

        assert!(!verifier.verify_eval_batched(&xs, &commitment, &witness))
    }

    #[test]
    fn test_eval_batched_all_points() {
        let params = test_setup::<15>();
        let (mut prover, verifier) = test_participants(&params);
        let polynomial = random_polynomial(8, 15);
        let commitment = prover.commit(polynomial.clone());

        let mut xs: Vec<Scalar> = Vec::with_capacity(polynomial.num_coeffs());
        let mut ys: Vec<Scalar> = Vec::with_capacity(polynomial.num_coeffs());
        for _ in 0..polynomial.num_coeffs() {
            let x: Scalar = RNG_1.lock().unwrap().gen::<u64>().into();
            xs.push(x);
            ys.push(polynomial.eval(x));
        }

        let witness = prover.create_witness_batched(xs.as_slice(), ys.as_slice()).unwrap();
        assert!(verifier.verify_eval_batched(xs.as_slice(), &commitment, &witness));
    }
}
