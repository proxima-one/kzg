use std::fmt::Debug;
use pairing::{
    group::{ff::Field, prime::PrimeCurveAffine, Curve},
    Engine,
};
use thiserror::Error;

#[cfg(feature = "serde_support")]
use serde::{Serialize, Deserialize};

pub mod wrapper_types;
pub mod utils;
pub mod ft;
pub mod polynomial;
pub mod worker;

use polynomial::{Polynomial, SubProductTree, op_tree};
use wrapper_types::{G1Affine, G2Affine};

/// parameters from tested setup
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde_support", derive(Serialize, Deserialize))]
pub struct KZGParams<E: Engine> {
    /// generator of g
    g: G1Affine<E>,
    /// generator of G2
    h: G2Affine<E>,
    /// g^alpha^1, g^alpha^2, ...
    gs: Vec<G1Affine<E>>,
    /// g^alpha^1, g^alpha^2, ...
    hs: Vec<G2Affine<E>>,
}

/// the commitment - "C" in the paper. It's a single group element
pub type KZGCommitment<E> = G1Affine<E>;
/// A witness for a single element - "w_i" in the paper. It's a group element.
pub type KZGWitness<E> = G1Affine<E>;

// A witness for a several elements - "w_B" in the paper. It's a single group element plus a polynomial
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KZGBatchWitness<E: Engine> {
    r: Polynomial<E::Fr>,
    w: G1Affine<E>,
}

impl<E: Engine> KZGBatchWitness<E> {
    pub fn elem(&self) -> G1Affine<E> {
        self.w
    }

    pub fn elem_ref(&self) -> &G1Affine<E> {
        &self.w
    }

    pub fn polynomial(&self) -> &Polynomial<E::Fr> {
        &self.r
    }

    pub fn from_inner(r: Polynomial<E::Fr>, w: G1Affine<E>) -> Self {
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
pub struct KZGProver<'params, E: Engine> {
    parameters: &'params KZGParams<E>,
    polynomial: Option<Polynomial<E::Fr>>,
    commitment: Option<KZGCommitment<E>>,
}

#[derive(Debug, Clone)]
pub struct KZGVerifier<'params, E: Engine> {
    parameters: &'params KZGParams<E>,
}

impl<'params, E: Engine> KZGProver<'params, E> {
    /// initializes `polynomial` to zero polynomial
    pub fn new(parameters: &'params KZGParams<E>) -> Self {
        Self {
            parameters,
            polynomial: None,
            commitment: None,
        }
    }

    pub fn parameters(&self) -> &'params KZGParams<E> {
        self.parameters
    }

    pub fn commit(&mut self, polynomial: Polynomial<E::Fr>) -> KZGCommitment<E> {
        let mut commitment = self.parameters.g.into_inner() * polynomial.coeffs[0];
        for i in 0..polynomial.degree() {
            commitment += self.parameters.gs[i].into_inner() * polynomial.coeffs[i + 1];
        }

        self.polynomial = Some(polynomial);
        let commitment = G1Affine::from_inner(commitment.to_affine());
        self.commitment = Some(commitment);
        commitment
    }

    pub fn open(&self) -> Result<Polynomial<E::Fr>, KZGError> {
        self.polynomial.clone().ok_or(KZGError::NoPolynomial)
    }

    pub fn commitment(&self) -> Option<KZGCommitment<E>> {
        self.commitment
    }

    pub fn commitment_ref(&self) -> Option<&KZGCommitment<E>> {
        self.commitment.as_ref()
    }

    pub fn set_commitment(&mut self, commitment: KZGCommitment<E>, polynomial: Polynomial<E::Fr>) {
        self.commitment = Some(commitment);
        self.polynomial = Some(polynomial);
    }

    pub fn has_commitment(&self) -> bool {
        self.commitment.is_some()
    }

    pub fn polynomial(&self) -> Option<&Polynomial<E::Fr>> {
        self.polynomial.as_ref()
    }

    pub fn has_polynomial(&self) -> bool {
        self.polynomial.is_some()
    }

    pub fn create_witness(&self, (x, y): (E::Fr, E::Fr)) -> Result<KZGWitness<E>, KZGError> {
        match self.polynomial {
            None => Err(KZGError::NoPolynomial),
            Some(ref polynomial) => {
                let mut dividend = polynomial.clone();
                dividend.coeffs[0] -= y;

                let divisor = Polynomial::new_from_coeffs(vec![-x, E::Fr::one()], 1);
                match dividend.long_division(&divisor) {
                    // by polynomial remainder theorem, if (x - point.x) does not divide self.polynomial, then
                    // self.polynomial(point.y) != point.1
                    (_, Some(_)) => Err(KZGError::PointNotOnPolynomial),
                    (psi, None) => {
                        let mut w = self.parameters.g.into_inner() * psi.coeffs[0];
                        for i in 0..psi.degree() {
                            w += self.parameters.gs[i].into_inner() * psi.coeffs[i + 1];
                        }

                        Ok(G1Affine::from_inner(w.to_affine()))
                    }
                }
            }
        }
    }

    pub fn create_witness_batched(
        &self,
        xs: &[E::Fr],
        ys: &[E::Fr]
    ) -> Result<KZGBatchWitness<E>, KZGError> {
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
                        let mut w = self.parameters.g.into_inner() * psi.coeffs[0];
                        for i in 0..psi.degree() {
                            w += self.parameters.gs[i].into_inner() * psi.coeffs[i + 1];
                        }
                        Ok(KZGBatchWitness {
                            r: interpolation,
                            w: G1Affine::from_inner(w.to_affine()),
                        })
                    }
                }
            }
        }
    }
}

impl<'params, E: Engine> KZGVerifier<'params, E> {
    pub fn new(parameters: &'params KZGParams<E>) -> Self {
        KZGVerifier { parameters }
    }

    pub fn verify_poly(
        &self,
        commitment: &KZGCommitment<E>,
        polynomial: &Polynomial<E::Fr>,
    ) -> bool {
        let mut check = self.parameters.g.into_inner() * polynomial.coeffs[0];
        for i in 0..polynomial.degree() {
            check += self.parameters.gs[i].into_inner() * polynomial.coeffs[i + 1];
        }

        check.to_affine() == commitment.into_inner()
    }

    pub fn verify_eval(
        &self,
        (x, y): (E::Fr, E::Fr),
        commitment: &KZGCommitment<E>,
        witness: &KZGWitness<E>,
    ) -> bool {
        let lhs = E::pairing(
            witness,
            &(self.parameters.hs[0].to_curve() - self.parameters.h.into_inner() * x).to_affine(),
        );
        let rhs = E::pairing(
            &(commitment.to_curve() - self.parameters.g.into_inner() * y).to_affine(),
            &self.parameters.h,
        );

        lhs == rhs
    }

    pub fn verify_eval_batched(
        &self,
        xs: &[E::Fr],
        commitment: &KZGCommitment<E>,
        witness: &KZGBatchWitness<E>,
    ) -> bool {
        let z: Polynomial<E::Fr> = op_tree(
            xs.len(),
            &|i| {
                let mut coeffs = vec![-xs[i], E::Fr::one()];
                coeffs[0] = -xs[i];
                coeffs[1] = E::Fr::one();
                Polynomial::new_from_coeffs(coeffs, 1)
            },
            &|a, b| a.best_mul(&b),
        );

        let mut hz = self.parameters.h.into_inner() * z.coeffs[0];
        for i in 0..z.degree() {
            hz += self.parameters.hs[i].into_inner() * z.coeffs[i + 1];
        }

        let mut gr = self.parameters.g.into_inner() * witness.r.coeffs[0];
        for i in 0..witness.r.degree() {
            gr += self.parameters.gs[i].into_inner() * witness.r.coeffs[i + 1];
        }

        let lhs = E::pairing(&witness.w, &hz.to_affine());
        let rhs = E::pairing(
            &(commitment.to_curve() - gr).to_affine(),
            &self.parameters.h,
        );

        lhs == rhs
    }
}

pub fn setup<E: Engine>(s: E::Fr, num_coeffs: usize) -> KZGParams<E> {
    let g = E::G1Affine::generator();
    let h = E::G2Affine::generator();

    let mut gs = vec![G1Affine::from_inner(g); num_coeffs];
    let mut hs = vec![G2Affine::from_inner(h); num_coeffs];

    let mut curr = g;

    for g in gs.iter_mut() {
        let inner = (curr * s).to_affine();
        *g = G1Affine::from_inner(inner);
        curr = inner;
    }

    let mut curr = h;
    for h in hs.iter_mut() {
        let inner = (curr * s).to_affine();
        *h = G2Affine::from_inner(inner);
        curr = inner;
    }

    KZGParams { g: G1Affine::from_inner(g), h: G2Affine::from_inner(h), gs, hs }
}

#[cfg(any(csprng_setup, test))]
use rand::random;

#[cfg(any(csprng_setup, test))]
pub fn csprng_setup<E: Engine>(num_coeffs: usize) -> KZGParams<E> {
    let s: E::Fr = random::<u64>().into();
    setup(s, num_coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use bls12_381::{Bls12, Scalar};
    use lazy_static::lazy_static;
    use pairing::group::ff::PrimeFieldBits;
    use rand::{rngs::SmallRng, Rng, SeedableRng};
    use std::sync::Mutex;

    const RNG_SEED_0: [u8; 32] = [42; 32];
    const RNG_SEED_1: [u8; 32] = [79; 32];

    lazy_static! {
        static ref RNG_0: Mutex<SmallRng> = Mutex::new(SmallRng::from_seed(RNG_SEED_0));
        static ref RNG_1: Mutex<SmallRng> = Mutex::new(SmallRng::from_seed(RNG_SEED_1));
    }

    fn test_setup<E: Engine, const MAX_COEFFS: usize>() -> KZGParams<E> {
        let s: E::Fr = RNG_0.lock().unwrap().gen::<u64>().into();
        setup(s, MAX_COEFFS)
    }

    fn test_participants<'params, E: Engine>(
        params: &'params KZGParams<E>,
    ) -> (KZGProver<'params, E>, KZGVerifier<'params, E>) {
        let prover = KZGProver::new(params);
        let verifier = KZGVerifier::new(params);

        (prover, verifier)
    }

    // never returns zero polynomial
    fn random_polynomial<S: PrimeFieldBits>(min_coeffs: usize, max_coeffs: usize) -> Polynomial<S> {
        let num_coeffs = RNG_1.lock().unwrap().gen_range(min_coeffs..max_coeffs);
        let mut coeffs = vec![S::zero(); max_coeffs];

        for i in 0..num_coeffs {
            coeffs[i] = RNG_1.lock().unwrap().gen::<u64>().into();
        }

        let mut poly = Polynomial::new_from_coeffs(coeffs, num_coeffs - 1);
        poly.shrink_degree();
        poly
    }

    fn assert_verify_poly<E: Engine + Debug>(
        verifier: &KZGVerifier<E>,
        commitment: &KZGCommitment<E>,
        polynomial: &Polynomial<E::Fr>,
    ) {
        assert!(
            verifier.verify_poly(&commitment, &polynomial),
            "verify_poly failed for commitment {:#?} and polynomial {:#?}",
            commitment,
            polynomial
        );
    }

    fn assert_verify_poly_fails<E: Engine + Debug>(
        verifier: &KZGVerifier<E>,
        commitment: &KZGCommitment<E>,
        polynomial: &Polynomial<E::Fr>,
    ) {
        assert!(
            !verifier.verify_poly(&commitment, &polynomial),
            "expected verify_poly to fail for commitment {:#?} and polynomial {:#?} but it didn't",
            commitment,
            polynomial
        );
    }

    fn assert_verify_eval<E: Engine + Debug>(
        verifier: &KZGVerifier<E>,
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

    fn assert_verify_eval_fails<E: Engine + Debug>(
        verifier: &KZGVerifier<E>,
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

        let polynomial = random_polynomial(1, 10);
        let commitment = prover.commit(polynomial.clone());

        assert_verify_poly(&verifier, &commitment, &polynomial);
        assert_verify_poly_fails(&verifier, &commitment, &random_polynomial(1, 10));
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

        let polynomial = random_polynomial(4, 8);
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

        let polynomial = random_polynomial(5, 13);
        let commitment = prover.commit(polynomial.clone());

        let x: Scalar = RNG_1.lock().unwrap().gen::<u64>().into();
        let y = polynomial.eval(x);

        let witness = prover.create_witness((x, y)).unwrap();
        assert_verify_eval(&verifier, (x, y), &commitment, &witness);

        let y_prime = random_field_elem_neq::<Bls12>(y);
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
        let params = test_setup::<Bls12, 15>();
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
        let params = test_setup::<Bls12, 15>();
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
