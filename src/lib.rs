use std::fmt::Debug;
use blstrs::{G1Affine, G1Projective, G2Affine, G2Projective, Scalar, pairing};
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
    /// g, g^alpha^1, g^alpha^2, ...
    gs: Vec<G1Affine>,
    /// h, h^alpha^1, h^alpha^2, ...
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
        let gs: Vec<G1Projective> = self.parameters.gs.iter().take(polynomial.num_coeffs()).map(|g| g.to_curve()).collect();
        let commitment = G1Projective::multi_exp(gs.as_slice(), polynomial.slice_coeffs());

        // let mut commitment = self.parameters.gs[0] * polynomial.coeffs[0];
        // for i in 1..polynomial.num_coeffs() {
        //     commitment += self.parameters.gs[i] * polynomial.coeffs[i];
        // }

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
                        let w = if psi.num_coeffs() == 1 {
                            self.parameters.gs[0] * psi.coeffs[0]
                        } else {
                            let gs: Vec<G1Projective> = self.parameters.gs.iter().take(psi.num_coeffs()).map(|g| g.to_curve()).collect();
                            G1Projective::multi_exp(gs.as_slice(), psi.slice_coeffs())
                        };

                        // let mut w = self.parameters.gs[0] * psi.coeffs[0];
                        // for i in 1..psi.num_coeffs() {
                        //     w += self.parameters.gs[i] * psi.coeffs[i];
                        // }

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
                        let w = if psi.num_coeffs() == 1 {
                            self.parameters.gs[0] * psi.coeffs[0]
                        } else {
                            let gs: Vec<G1Projective> = self.parameters.gs.iter().take(psi.num_coeffs()).map(|g| g.to_curve()).collect();
                            G1Projective::multi_exp(gs.as_slice(), psi.slice_coeffs())
                        };

                        // let mut w = self.parameters.gs[0] * psi.coeffs[0];
                        // for i in 1..psi.num_coeffs() {
                        //     w += self.parameters.gs[i] * psi.coeffs[i];
                        // }

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
        let gs: Vec<G1Projective> = self.parameters.gs.iter().take(polynomial.num_coeffs()).map(|g| g.to_curve()).collect();
        let check = G1Projective::multi_exp(gs.as_slice(), polynomial.slice_coeffs());

        // let mut check = self.parameters.gs[0] * polynomial.coeffs[0];
        // for i in 1..polynomial.num_coeffs() {
        //     check += self.parameters.gs[i] * polynomial.coeffs[i];
        // }

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
            &(self.parameters.hs[1].to_curve() - self.parameters.hs[0] * x).to_affine(),
        );
        let rhs = pairing(
            &(commitment.to_curve() - self.parameters.gs[0] * y).to_affine(),
            &self.parameters.hs[0],
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

        // let mut hz = self.parameters.hs[0] * z.coeffs[0];
        // for i in 1..z.num_coeffs() {
        //     hz += self.parameters.hs[i] * z.coeffs[i];
        // }

        let hz = if z.num_coeffs() == 1 {
            self.parameters.hs[0] * z.coeffs[0]
        } else {
            let hs: Vec<G2Projective> = self.parameters.hs.iter().take(z.num_coeffs()).map(|h| h.to_curve()).collect();
            G2Projective::multi_exp(hs.as_slice(), z.slice_coeffs())
        };


        // let mut gr= self.parameters.gs[0] * witness.r.coeffs[0];
        // for i in 1..witness.r.num_coeffs() {
        //     gr += self.parameters.gs[i] * witness.r.coeffs[i];
        // }

        let gr = if witness.r.num_coeffs() == 1 {
            self.parameters.gs[0] * witness.r.coeffs[0]
        } else {
            let gs: Vec<G1Projective> = self.parameters.gs.iter().take(witness.r.num_coeffs()).map(|g| g.to_curve()).collect();
            G1Projective::multi_exp(gs.as_slice(), witness.r.slice_coeffs())
        };

        let lhs = pairing(&witness.w, &hz.to_affine());
        let rhs = pairing(
            &(commitment.to_curve() - gr).to_affine(),
            &self.parameters.hs[0],
        );

        lhs == rhs
    }
}

pub fn setup(s: Scalar, num_coeffs: usize) -> KZGParams {
    let mut gs = vec![G1Affine::generator(); num_coeffs + 1];
    let mut hs = vec![G2Affine::generator(); num_coeffs + 1];

    let mut curr = gs[0];
    for g in gs.iter_mut().skip(1) {
        *g = (curr * s).to_affine();
        curr = *g;
    }

    let mut curr = hs[0];
    for h in hs.iter_mut().skip(1) {
        *h = (curr * s).to_affine();
        curr = *h;
    }

    KZGParams { gs, hs }
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
    use rand::{rngs::SmallRng, Rng, SeedableRng};

    const RNG_SEED: [u8; 32] = [69; 32];

    fn test_setup<const MAX_COEFFS: usize>(rng: &mut SmallRng) -> KZGParams {
        let s: Scalar = rng.gen::<u64>().into();
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
    fn random_polynomial(rng: &mut SmallRng, min_coeffs: usize, max_coeffs: usize) -> Polynomial {
        let num_coeffs = rng.gen_range(min_coeffs..max_coeffs);
        let mut coeffs = vec![Scalar::zero(); max_coeffs];

        for i in 0..num_coeffs {
            coeffs[i] = rng.gen::<u64>().into();
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
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<10>(&mut rng);

        let (mut prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(&mut rng, 2, 12);
        let commitment = prover.commit(polynomial.clone());

        assert_verify_poly(&verifier, &commitment, &polynomial);
        assert_verify_poly_fails(&verifier, &commitment, &random_polynomial(&mut rng, 2, 12));
    }

    fn random_field_elem_neq(val: Scalar) -> Scalar {
        let mut rng = SmallRng::from_seed(RNG_SEED); 
        let mut v: Scalar = rng.gen::<u64>().into();
        while v == val {
            v = rng.gen::<u64>().into();
        }

        v
    }

    #[test]
    fn test_modify_single_coeff() {
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<8>(&mut rng);

        let (mut prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(&mut rng, 3, 8);
        let commitment = prover.commit(polynomial.clone());

        let mut modified_polynomial = polynomial.clone();
        let new_coeff = random_field_elem_neq(modified_polynomial.coeffs[2]);
        modified_polynomial.coeffs[2] = new_coeff;

        assert_verify_poly(&verifier, &commitment, &polynomial);
        assert_verify_poly_fails(&verifier, &commitment, &modified_polynomial);
    }

    #[test]
    fn test_eval_basic() {
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<13>(&mut rng);

        let (mut prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(&mut rng, 5, 13);
        let commitment = prover.commit(polynomial.clone());

        let x: Scalar = rng.gen::<u64>().into();
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
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<15>(&mut rng);

        let (mut prover, verifier) = test_participants(&params);
        let polynomial = random_polynomial(&mut rng, 8, 15);
        let commitment = prover.commit(polynomial.clone());

        let mut xs: Vec<Scalar> = Vec::with_capacity(8);
        let mut ys: Vec<Scalar> = Vec::with_capacity(8);
        for _ in 0..8 {
            let x: Scalar = rng.gen::<u64>().into();
            xs.push(x);
            ys.push(polynomial.eval(x));
        }

        let witness = prover.create_witness_batched(xs.as_slice(), ys.as_slice()).unwrap();
        assert!(verifier.verify_eval_batched(xs.as_slice(), &commitment, &witness));

        let mut xs: Vec<Scalar> = Vec::with_capacity(8);
        let mut ys: Vec<Scalar> = Vec::with_capacity(8);
        for _ in 0..8 {
            let x: Scalar = rng.gen::<u64>().into();
            xs.push(x);
            ys.push(polynomial.eval(x));
        }

        assert!(!verifier.verify_eval_batched(&xs, &commitment, &witness))
    }

    #[test]
    fn test_eval_batched_all_points() {
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<15>(&mut rng);
        
        let (mut prover, verifier) = test_participants(&params);
        let polynomial = random_polynomial(&mut rng, 13, 14);
        let commitment = prover.commit(polynomial.clone());
        
        let mut xs: Vec<Scalar> = Vec::with_capacity(polynomial.num_coeffs());
        let mut ys: Vec<Scalar> = Vec::with_capacity(polynomial.num_coeffs());
        for _ in 0..polynomial.num_coeffs() {
            let x: Scalar = rng.gen::<u64>().into();
            xs.push(x);
            ys.push(polynomial.eval(x));
        }

        let witness = prover.create_witness_batched(xs.as_slice(), ys.as_slice()).unwrap();
        assert!(verifier.verify_eval_batched(xs.as_slice(), &commitment, &witness));
    }
}
