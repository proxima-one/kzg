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
pub struct KZGBatchWitness<const N: usize> {
    r: Polynomial<N>,
    w: G1Affine,
}

impl<const N: usize> KZGBatchWitness<N> {
    pub fn elem(&self) -> G1Affine {
        self.w
    }

    pub fn elem_ref(&self) -> &G1Affine {
        &self.w
    }

    pub fn polynomial(&self) -> &Polynomial<N> {
        &self.r
    }

    pub fn new(r: Polynomial<N>, w: G1Affine) -> Self {
        KZGBatchWitness { r, w }
    }
}

#[derive(Debug, Clone)]
pub struct KZGProver<'params, const N: usize> {
    parameters: &'params KZGParams<N>,
}

#[derive(Debug, Clone)]
pub struct KZGVerifier<'params, const N: usize> {
    parameters: &'params KZGParams<N>,
}

impl<'params, const N: usize> KZGProver<'params, N> {
    /// initializes `polynomial` to zero polynomial
    pub fn new(parameters: &'params KZGParams<N>) -> Self {
        Self {
            parameters,
        }
    }

    pub fn parameters(&self) -> &'params KZGParams<N> {
        self.parameters
    }

    pub fn commit(&self, polynomial: &Polynomial<N>) -> KZGCommitment {
        let gs = &self.parameters.gs[..polynomial.num_coeffs()];
        let commitment = G1Projective::multi_exp(gs, polynomial.slice_coeffs());

        commitment.to_affine()
    }

    pub fn create_witness(&self, polynomial: &Polynomial<N>, (x, y): (Scalar, Scalar)) -> Result<KZGWitness, KZGError> {
        let mut dividend = polynomial.clone();
        dividend.coeffs[0] -= y;

        let mut divisor_coeffs = [Scalar::zero(); N];
        divisor_coeffs[0] = -x;
        divisor_coeffs[1] = Scalar::one();
        let divisor = Polynomial::new_from_coeffs(divisor_coeffs, 1);

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
        polynomial: &Polynomial<N>,
        xs: &[Scalar],
        ys: &[Scalar],
    ) -> Result<KZGBatchWitness<N>, KZGError> {
        assert!(xs.len() == ys.len());
        assert!(xs.len() <= polynomial.num_coeffs());

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

impl<'params, const N: usize> KZGVerifier<'params, N> {
    pub fn new(parameters: &'params KZGParams<N>) -> Self {
        KZGVerifier { parameters }
    }

    pub fn verify_poly(&self, commitment: &KZGCommitment, polynomial: &Polynomial<N>) -> bool {
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
        witness: &KZGBatchWitness<N>,
    ) -> bool {
        let z: Polynomial<N> = op_tree(
            xs.len(),
            &|i| {
                let mut coeffs = [Scalar::zero(); N];
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

    fn test_setup<const N: usize>(rng: &mut SmallRng) -> KZGParams<N> {
        let s: Scalar = rng.gen::<u64>().into();
        setup::<N>(s)
    }

    fn test_participants<'params, const N: usize>(
        params: &'params KZGParams<N>,
    ) -> (KZGProver<'params, N>, KZGVerifier<'params, N>) {
        let prover = KZGProver::new(params);
        let verifier = KZGVerifier::new(params);

        (prover, verifier)
    }

    // never returns zero polynomial
    fn random_polynomial<const N: usize>(rng: &mut SmallRng, min_coeffs: usize) -> Polynomial<N> {
        let num_coeffs = rng.gen_range(min_coeffs..N);
        let mut coeffs = [Scalar::zero(); N];

        for i in 0..num_coeffs {
            coeffs[i] = rng.gen::<u64>().into();
        }

        let mut poly = Polynomial::new_from_coeffs(coeffs, num_coeffs - 1);
        poly.shrink_degree();
        poly
    }

    fn assert_verify_poly<const N: usize>(
        verifier: &KZGVerifier<N>,
        commitment: &KZGCommitment,
        polynomial: &Polynomial<N>,
    ) {
        assert!(
            verifier.verify_poly(&commitment, &polynomial),
            "verify_poly failed for commitment {:#?} and polynomial {:#?}",
            commitment,
            polynomial
        );
    }

    fn assert_verify_poly_fails<const N: usize>(
        verifier: &KZGVerifier<N>,
        commitment: &KZGCommitment,
        polynomial: &Polynomial<N>,
    ) {
        assert!(
            !verifier.verify_poly(&commitment, &polynomial),
            "expected verify_poly to fail for commitment {:#?} and polynomial {:#?} but it didn't",
            commitment,
            polynomial
        );
    }

    fn assert_verify_eval<const N: usize>(
        verifier: &KZGVerifier<N>,
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

    fn assert_verify_eval_fails<const N: usize>(
        verifier: &KZGVerifier<N>,
        point: (Scalar, Scalar),
        commitment: &KZGCommitment,
        witness: &KZGWitness,
    ) {
        assert!(!verifier.verify_eval(point, &commitment, &witness), "expected verify_eval to fail for for point {:#?}, commitment {:#?}, and witness {:#?}, but it didn't", point, commitment, witness);
    }

    #[test]
    fn test_basic() {
        const N: usize = 16;

        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<N>(&mut rng);

        let (prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial::<N>(&mut rng, 2);
        let commitment = prover.commit(&polynomial);

        assert_verify_poly(&verifier, &commitment, &polynomial);
        assert_verify_poly_fails(&verifier, &commitment, &random_polynomial(&mut rng, 2));
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
        const N: usize = 8;
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<N>(&mut rng);

        let (prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(&mut rng, 3);
        let commitment = prover.commit(&polynomial);

        let mut modified_polynomial = polynomial.clone();
        let new_coeff = random_field_elem_neq(modified_polynomial.coeffs[2]);
        modified_polynomial.coeffs[2] = new_coeff;

        assert_verify_poly(&verifier, &commitment, &polynomial);
        assert_verify_poly_fails(&verifier, &commitment, &modified_polynomial);
    }

    #[test]
    fn test_eval_basic() {
        const N: usize = 16;
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<N>(&mut rng);

        let (prover, verifier) = test_participants(&params);

        let polynomial = random_polynomial(&mut rng, 5);
        let commitment = prover.commit(&polynomial);

        let x: Scalar = rng.gen::<u64>().into();
        let y = polynomial.eval(x);

        let witness = prover.create_witness(&polynomial, (x, y)).unwrap();
        assert_verify_eval(&verifier, (x, y), &commitment, &witness);

        let y_prime = random_field_elem_neq(y);
        assert_verify_eval_fails(&verifier, (x, y_prime), &commitment, &witness);

        // test degree 1 edge case
        let mut coeffs = [Scalar::zero(); N];
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
        const N: usize = 8;
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<N>(&mut rng);

        let (prover, verifier) = test_participants(&params);
        let polynomial = random_polynomial(&mut rng, 5);
        let commitment = prover.commit(&polynomial);

        let mut xs = [Scalar::zero(); N];
        let mut ys  = [Scalar::zero(); N];
        for i in 0..N {
            let x: Scalar = rng.gen::<u64>().into();
            xs[i] = x;
            ys[i] = polynomial.eval(x);
        }

        let witness = prover
            .create_witness_batched(&polynomial, &xs[..5], &ys[..5])
            .unwrap();
        assert!(verifier.verify_eval_batched(&xs, &commitment, &witness));

        let mut xs = [Scalar::zero(); N];
        let mut ys = [Scalar::zero(); N];
        for i in 0..N {
            let x: Scalar = rng.gen::<u64>().into();
            xs[i] = x;
            ys[i] = polynomial.eval(x);
        }

        assert!(!verifier.verify_eval_batched(&xs[..5], &commitment, &witness))
    }

    #[test]
    fn test_eval_batched_all_points() {
        const N: usize = 16;
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<N>(&mut rng);

        let (prover, verifier) = test_participants(&params);
        let polynomial = random_polynomial(&mut rng, 13);
        let commitment = prover.commit(&polynomial);

        let mut xs = [Scalar::zero(); N];
        let mut ys = [Scalar::zero(); N];
        for i in 0..polynomial.num_coeffs() {
            let x: Scalar = rng.gen::<u64>().into();
            xs[i] = x;
            ys[i] = polynomial.eval(x);
        }

        let witness = prover
            .create_witness_batched(&polynomial, &xs[..], &ys[..])
            .unwrap();
        assert!(verifier.verify_eval_batched(&xs[..], &commitment, &witness));
    }
}
