use blstrs::{pairing, G1Affine, G1Projective, G2Projective, Scalar};
use pairing::group::{ff::Field, prime::PrimeCurveAffine, Curve};
use std::fmt::Debug;

#[cfg(feature = "serde_support")]
use serde::{Deserialize, Serialize};

use crate::polynomial::{op_tree, Polynomial, SubProductTree};
use crate::{KZGCommitment, KZGError, KZGParams, KZGWitness};

// A witness for a several elements - "w_B" in the paper. It's a single group element plus a polynomial
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde_support", derive(Serialize, Deserialize))]
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

    pub fn new(r: Polynomial, w: G1Affine) -> Self {
        KZGBatchWitness { r, w }
    }
}

#[derive(Debug, Clone)]
pub struct KZGProver<'params> {
    parameters: &'params KZGParams,
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
        }
    }

    pub fn parameters(&self) -> &'params KZGParams {
        self.parameters
    }

    pub fn commit(&self, polynomial: &Polynomial) -> KZGCommitment {
        let gs = &self.parameters.gs[..polynomial.num_coeffs()];
        let commitment = G1Projective::multi_exp(gs, polynomial.slice_coeffs());

        commitment.to_affine()
    }

    pub fn create_witness(&self, polynomial: &Polynomial, (x, y): (Scalar, Scalar)) -> Result<KZGWitness, KZGError> {
        let mut dividend = polynomial.clone();
        dividend.coeffs[0] -= y;

        let divisor = Polynomial::new_from_coeffs(vec![-x, Scalar::one()], 1);
        match dividend.long_division(&divisor) {
            // by polynomial remainder theorem, if (x - point.x) does not divide self.polynomial, then
            // self.polynomial(point.y) != point.1
            (_, Some(_)) => Err(KZGError::PointNotOnPolynomial),
            (psi, None) if psi.num_coeffs() == 1 => Ok((self.parameters.gs[0] * psi.coeffs[0]).to_affine()),
            (psi, None) => {
                let gs = &self.parameters.gs[..psi.num_coeffs()];
                Ok(G1Projective::multi_exp(gs, psi.slice_coeffs()).to_affine())
            }
        }
    }

    pub fn create_witness_batched(
        &self,
        polynomial: &Polynomial,
        xs: &[Scalar],
        ys: &[Scalar],
    ) -> Result<KZGBatchWitness, KZGError> {
        let tree = SubProductTree::new_from_points(xs);

        let interpolation = Polynomial::lagrange_interpolation_with_tree(xs, ys, &tree);

        let numerator = polynomial - &interpolation;
        let (psi, rem) = numerator.long_division(&tree.product);
        match rem {
            Some(_) => Err(KZGError::PointNotOnPolynomial),
            None => {
                let w = if psi.num_coeffs() == 1 {
                    self.parameters.gs[0] * psi.coeffs[0]
                } else {
                    let gs = &self.parameters.gs[..psi.num_coeffs()];
                    G1Projective::multi_exp(gs, psi.slice_coeffs())
                };

                Ok(KZGBatchWitness {
                    r: interpolation,
                    w: w.to_affine(),
                })
            }
        }
    }
}

impl<'params> KZGVerifier<'params> {
    pub fn new(parameters: &'params KZGParams) -> Self {
        KZGVerifier { parameters }
    }

    pub fn verify_poly(&self, commitment: &KZGCommitment, polynomial: &Polynomial) -> bool {
        let gs = &self.parameters.gs[..polynomial.num_coeffs()];
        let check = G1Projective::multi_exp(gs, polynomial.slice_coeffs());

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
            &(self.parameters.hs[1] - self.parameters.hs[0] * x).to_affine(),
        );
        let rhs = pairing(
            &(commitment.to_curve() - self.parameters.gs[0] * y).to_affine(),
            &self.parameters.hs[0].to_affine(),
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

        let hz = if z.num_coeffs() == 1 {
            self.parameters.hs[0] * z.coeffs[0]
        } else {
            let hs = &self.parameters.hs[..z.num_coeffs()];
            G2Projective::multi_exp(hs, z.slice_coeffs())
        };

        let gr = if witness.r.num_coeffs() == 1 {
            self.parameters.gs[0] * witness.r.coeffs[0]
        } else {
            let gs = &self.parameters.gs[..witness.r.num_coeffs()];
            G1Projective::multi_exp(gs, witness.r.slice_coeffs())
        };

        let lhs = pairing(&witness.w, &hz.to_affine());
        let rhs = pairing(
            &(commitment.to_curve() - gr).to_affine(),
            &self.parameters.hs[0].to_affine(),
        );

        lhs == rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::setup;
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
        let params = test_setup::<12>(&mut rng);

        let (prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(&mut rng, 2, 12);
        let commitment = prover.commit(&polynomial);

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

        let (prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(&mut rng, 3, 8);
        let commitment = prover.commit(&polynomial);

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

        let (prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(&mut rng, 5, 13);
        let commitment = prover.commit(&polynomial);

        let x: Scalar = rng.gen::<u64>().into();
        let y = polynomial.eval(x);

        let witness = prover.create_witness(&polynomial, (x, y)).unwrap();
        assert_verify_eval(&verifier, (x, y), &commitment, &witness);

        let y_prime = random_field_elem_neq(y);
        assert_verify_eval_fails(&verifier, (x, y_prime), &commitment, &witness);

        // test degree 1 edge case
        let mut coeffs = vec![Scalar::zero(); 13];
        coeffs[0] = 3.into();
        coeffs[1] = 1.into();
        let polynomial = Polynomial::new(coeffs);

        let commitment = prover.commit(&polynomial);
        let witness = prover.create_witness(&polynomial, (1.into(), 4.into())).unwrap();
        assert_verify_eval(&verifier, (1.into(), 4.into()), &commitment, &witness);
        assert_verify_eval_fails(&verifier, (1.into(), 5.into()), &commitment, &witness);
    }

    #[test]
    fn test_eval_batched() {
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<15>(&mut rng);

        let (prover, verifier) = test_participants(&params);
        let polynomial = random_polynomial(&mut rng, 8, 15);
        let commitment = prover.commit(&polynomial);

        let mut xs: Vec<Scalar> = Vec::with_capacity(8);
        let mut ys: Vec<Scalar> = Vec::with_capacity(8);
        for _ in 0..8 {
            let x: Scalar = rng.gen::<u64>().into();
            xs.push(x);
            ys.push(polynomial.eval(x));
        }

        let witness = prover
            .create_witness_batched(&polynomial, xs.as_slice(), ys.as_slice())
            .unwrap();
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

        let (prover, verifier) = test_participants(&params);
        let polynomial = random_polynomial(&mut rng, 13, 14);
        let commitment = prover.commit(&polynomial);

        let mut xs: Vec<Scalar> = Vec::with_capacity(polynomial.num_coeffs());
        let mut ys: Vec<Scalar> = Vec::with_capacity(polynomial.num_coeffs());
        for _ in 0..polynomial.num_coeffs() {
            let x: Scalar = rng.gen::<u64>().into();
            xs.push(x);
            ys.push(polynomial.eval(x));
        }

        let witness = prover
            .create_witness_batched(&polynomial, xs.as_slice(), ys.as_slice())
            .unwrap();
        assert!(verifier.verify_eval_batched(xs.as_slice(), &commitment, &witness));
    }
}
