use blstrs::{pairing, G1Affine, G1Projective, G2Affine, G2Projective, Scalar};
use pairing::group::ff::PrimeField;
use pairing::group::{ff::Field, Group, prime::PrimeCurveAffine, Curve};
use std::fmt::Debug;
use arrayref::array_ref;

#[cfg(feature = "serde_support")]
use serde::{Deserialize, Serialize};

use crate::ft::{EvaluationDomain, compute_omega};
use crate::polynomial::Polynomial;
use crate::{KZGCommitment, KZGError, KZGParams, KZGWitness};

// A witness for a several elements - "w_B" in the paper. It's a single group element plus a polynomial
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde_support", derive(Serialize, Deserialize))]
pub struct KZGBatchWitnessEvalForm<const N: usize> {
    r: EvaluationDomain<N>,
    w: G1Affine,
}

impl<const N: usize> KZGBatchWitnessEvalForm<N> {
    pub fn elem(&self) -> G1Affine {
        self.w
    }

    pub fn elem_ref(&self) -> &G1Affine {
        &self.w
    }

    pub fn polynomial(&self) -> &EvaluationDomain<N> {
        &self.r
    }

    pub fn new(r: EvaluationDomain<N>, w: G1Affine) -> Self {
        KZGBatchWitnessEvalForm { r, w }
    }
}

#[derive(Debug, Clone)]
pub struct KZGProverEvalForm<'params, const N: usize> {
    parameters: &'params KZGParams<N>,
    lagrange_basis_g: &'params [G1Projective; N],
    lagrange_basis_h: &'params [G2Projective; N],
    exp: u32,
    omega: Scalar,
}

#[derive(Debug, Clone)]
pub struct KZGVerifierEvalForm<'params, const N: usize> {
    exp: u32,
    omega: Scalar,
    parameters: &'params KZGParams<N>,
    lagrange_basis_g: &'params [G1Projective; N],
    lagrange_basis_h: &'params [G2Projective; N],
}

fn div_by_omega_i<const N: usize>(evals: &EvaluationDomain<N>, m: usize) -> EvaluationDomain<N> {
    let mut coeffs = [Scalar::zero(); N];
    let omega_m = evals.omega.pow_vartime(&[m as u64]);
    for (j, &f) in evals.coeffs.iter().enumerate() {
        if j == m {
            let mut qm = Scalar::zero();
            let am = Scalar::from(N as u64) * omega_m.invert().unwrap();

            for i in (0..N).filter(|&i| i != m) {
                let omega_i = evals.omega.pow_vartime(&[i as u64]);
                let ai = Scalar::from(N as u64) * omega_i.invert().unwrap();

                let mut term = evals.coeffs[i];
                term *= am * ai.invert().unwrap();
                term *= (omega_m - omega_i).invert().unwrap();
                qm += term;
            }

            coeffs[j] = qm;
        } else {
            let omega_j = evals.omega.pow_vartime(&[j as u64]);
            coeffs[j] = f * (omega_j - omega_m).invert().unwrap();
        }
    }

    evals.clone_with_different_coeffs(coeffs)
}

impl<'params, const N: usize> KZGProverEvalForm<'params, N> {
    /// initializes `polynomial` to zero polynomial
    pub fn new(
        parameters: &'params KZGParams<N>,
        lagrange_basis_g: &'params [G1Projective; N],
        lagrange_basis_h: &'params [G2Projective; N],
    ) -> Self {
        let (_, exp, omega) = compute_omega(parameters.gs.len()).unwrap();
        Self {
            parameters,
            lagrange_basis_g,
            lagrange_basis_h,
            exp,
            omega,
        }
    }

    pub fn parameters(&self) -> &'params KZGParams<N> {
        self.parameters
    }

    pub fn omega(&self) -> Scalar {
        self.omega
    }

    pub fn commit(&self, evals: &EvaluationDomain<N>) -> KZGCommitment {
        let gs = &self.lagrange_basis_g[..evals.len()];
        let commitment = G1Projective::multi_exp(gs, evals.as_ref());

        let commitment = commitment.to_affine();
        commitment
    }

    pub fn create_witness(&self, evals: &EvaluationDomain<N>, i: usize) -> KZGWitness {
        let y = evals.coeffs[i];
        let mut numerator = evals.clone();
        for coeff in numerator.coeffs.iter_mut() {
            *coeff -= y;
        }

        let q = div_by_omega_i(&numerator, i);
        let w = if q.coeffs.len() == 1 {
            self.lagrange_basis_g[0] * q.coeffs[0]
        } else {
            let gs = &self.lagrange_basis_g[..q.len()];
            G1Projective::multi_exp(gs, q.as_ref())
        };

        w.to_affine()
    }

    pub fn create_witness_all(&self) -> KZGWitness {
        // this should get turned into a constant by the compiler
        let w: G1Projective = G1Projective::identity() * Scalar::zero();
        w.to_affine()
    }
}

impl<'params, const N: usize> KZGVerifierEvalForm<'params, N> {
    pub fn new(parameters: &'params KZGParams<N>, lagrange_basis_g: &'params [G1Projective; N], lagrange_basis_h: &'params [G2Projective; N]) -> Self {
        let (_, exp, omega) = compute_omega(parameters.gs.len()).unwrap();
        KZGVerifierEvalForm {
            parameters,
            exp,
            omega,
            lagrange_basis_g,
            lagrange_basis_h
        }
    }

    pub fn verify_poly(&self, commitment: &KZGCommitment, evals: &EvaluationDomain<N>) -> bool {
        let mut evals = evals.clone();
        evals.ifft();
        let polynomial: Polynomial<N> = evals.into();

        let gs = &self.parameters.gs[..polynomial.num_coeffs()];
        let check = G1Projective::multi_exp(gs, polynomial.slice_coeffs());

        check.to_affine() == *commitment
    }

    pub fn verify_eval(
        &self,
        (i, y): (usize, Scalar),
        commitment: &KZGCommitment,
        witness: &KZGWitness,
    ) -> bool {
        let omega = Scalar::root_of_unity().pow_vartime(&[1 << (Scalar::S - self.exp)]);
        let omega_i = omega.pow_vartime(&[i as u64]);
        let lhs = pairing(
            witness,
            &(self.parameters.hs[1] - self.parameters.hs[0] * omega_i).to_affine(),
        );
        let rhs = pairing(
            &(commitment.to_curve() - self.parameters.gs[0] * y).to_affine(),
            &self.parameters.hs[0].to_affine(),
        );

        lhs == rhs
    }

    pub fn verify_eval_all(
        &self,
        ys: &[Scalar; N],
        commitment: &KZGCommitment,
        witness: &KZGWitness,
    ) -> bool {
        let mut z = EvaluationDomain::new([Scalar::zero(); N], self.exp, self.omega);
        z.coeffs[0] = -Scalar::one();
        z.coeffs[N - 1] = Scalar::one();

        let hs = &self.lagrange_basis_h[..z.len()];
        let hz = G2Projective::multi_exp(hs, &z.coeffs);

        let r = EvaluationDomain::new(ys.clone(), self.exp, self.omega);

        let gs = &self.lagrange_basis_g[..r.len()];
        let gr = G1Projective::multi_exp(gs, &r.coeffs);

        let lhs = pairing(witness, &hz.to_affine());
        let rhs = pairing(
            &(commitment.to_curve() - gr).to_affine(),
            &self.parameters.hs[0].to_affine(),
        );

        lhs == rhs
    }
}

/// params.gs.len() must be a power of two.
pub fn compute_lagrange_basis<const N: usize>(params: &KZGParams<N>) -> ([G1Projective; N], [G2Projective; N]) {
    let d = params.gs.len();
    assert!(d & (d - 1) == 0);

    let (d, _, omega) = compute_omega(params.gs.len()).unwrap();
    let mut gs = [G1Projective::identity(); N];
    let mut hs = [G2Projective::identity(); N];

    for i in 0..d {
        let xi = omega.pow_vartime(&[i as u64]);
        let mut l = Polynomial::new_single_degree_zero(Scalar::one());
        for j in (0..d).filter(|&j| j != i) {
            let xj = omega.pow_vartime(&[j as u64]);
            let mut coeffs = [Scalar::zero(); N];
            coeffs[0] = -xj;
            coeffs[1] = Scalar::one();
            l = l.best_mul(&Polynomial::new_from_coeffs(coeffs, 1));
            l = l.scalar_multiplication((xi - xj).invert().unwrap());
        }

        let coeffs = l.coeffs();
        let g = &params.gs[..coeffs.len()];
        gs[i] = G1Projective::multi_exp(g, &coeffs);

        let h = &params.hs[..coeffs.len()];
        hs[i] = G2Projective::multi_exp(h, &coeffs);
    }

    (gs, hs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::setup;
    use rand::{rngs::SmallRng, Rng, SeedableRng};

    const RNG_SEED: [u8; 32] = [69; 32];

    fn test_setup<const N: usize>(rng: &mut SmallRng) -> KZGParams<N> {
        assert!(N & (N - 1) == 0);

        let s: Scalar = rng.gen::<u64>().into();
        setup::<N>(s)
    }

    fn test_participants<'params, const N: usize>(
        params: &'params KZGParams<N>,
        lagrange_basis_g: &'params [G1Projective; N],
        lagrange_basis_h: &'params [G2Projective; N],
    ) -> (KZGProverEvalForm<'params, N>, KZGVerifierEvalForm<'params, N>) {
        let prover = KZGProverEvalForm::new(params, lagrange_basis_g, lagrange_basis_h);
        let verifier = KZGVerifierEvalForm::new(params, lagrange_basis_g, lagrange_basis_h);

        (prover, verifier)
    }

    fn random_evals<const N: usize>(rng: &mut SmallRng) -> EvaluationDomain<N> {
        let mut coeffs = [Scalar::zero(); N];

        for i in 0..N {
            coeffs[i] = rng.gen::<u64>().into();
        }

        EvaluationDomain::from_coeffs(coeffs).unwrap()
    }

    #[test]
    fn test_div_by_omega_i() {
        const N: usize = 16;
        let (_, exp, omega) = compute_omega(N).unwrap();
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let omega_3 = omega.pow_vartime(&[3 as u64]);

        let mut top = Polynomial::new(random_evals::<N>(&mut rng).coeffs);
        let y = top.eval(omega_3);
        top.coeffs[0] -= y;

        let mut bottom_coeffs = [Scalar::zero(); N];
        bottom_coeffs[0] = -omega_3;
        bottom_coeffs[1] = Scalar::one();
        let bottom = Polynomial::new(bottom_coeffs);

        let (naive, rem) = top.long_division(&bottom);
        assert!(rem.is_none());

        let mut top = EvaluationDomain::new(top.coeffs, exp, omega);
        top.fft();
        let mut smart = div_by_omega_i(&top, 3);
        smart.ifft();
        let smart: Polynomial<N> = smart.into();

        assert_eq!(smart, naive);
    }

    #[test]
    fn test_subtract_by_scalar() {
        const N: usize = 16;
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let mut evals = random_evals::<N>(&mut rng);

        let subend = evals.clone_with_different_coeffs([2.into(); N]);

        let mut naive = evals.clone();
        naive.ifft();
        let mut naive: Polynomial<N> = naive.into();
        naive.coeffs[0] -= Scalar::from(2);

        evals.sub_assign(&subend);
        evals.ifft();
        let smart: Polynomial<N> = evals.into();
        assert_eq!(smart, naive);
    }

    fn assert_verify_poly<const N: usize>(
        verifier: &KZGVerifierEvalForm<N>,
        commitment: &KZGCommitment,
        evals: &EvaluationDomain<N>,
    ) {
        assert!(
            verifier.verify_poly(&commitment, &evals),
            "verify_poly failed for commitment {:#?} and polynomial {:#?}",
            commitment,
            evals
        );
    }

    fn assert_verify_poly_fails<const N: usize>(
        verifier: &KZGVerifierEvalForm<N>,
        commitment: &KZGCommitment,
        evals: &EvaluationDomain<N>,
    ) {
        assert!(
            !verifier.verify_poly(&commitment, &evals),
            "expected verify_poly to fail for commitment {:#?} and polynomial {:#?} but it didn't",
            commitment,
            evals
        );
    }

    #[test]
    fn test_basic() {
        const N: usize = 16;
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<N>(&mut rng);
        let lagrange_basis = compute_lagrange_basis(&params);

        let (mut prover, verifier) = test_participants(&params, &lagrange_basis.0, &lagrange_basis.1);

        let evals = random_evals::<N>(&mut rng);
        let commitment = prover.commit(&evals);

        assert_verify_poly(&verifier, &commitment, &evals);
        assert_verify_poly_fails(&verifier, &commitment, &random_evals::<N>(&mut rng));
    }

    fn random_field_elem_neq(rng: &mut SmallRng, val: Scalar) -> Scalar {
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
        let lagrange_basis = compute_lagrange_basis(&params);

        let (prover, verifier) = test_participants(&params, &lagrange_basis.0, &lagrange_basis.1);

        let evals = random_evals::<N>(&mut rng);
        let commitment = prover.commit(&evals);

        let mut modified_evals = evals.clone();
        let new_coeff = random_field_elem_neq(&mut rng, modified_evals.coeffs[2]);
        modified_evals.coeffs[2] = new_coeff;

        assert_verify_poly(&verifier, &commitment, &evals);
        assert_verify_poly_fails(&verifier, &commitment, &modified_evals);
    }

    fn assert_verify_eval<const N: usize>(
        verifier: &KZGVerifierEvalForm<N>,
        point: (usize, Scalar),
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
        verifier: &KZGVerifierEvalForm<N>,
        point: (usize, Scalar),
        commitment: &KZGCommitment,
        witness: &KZGWitness,
    ) {
        assert!(!verifier.verify_eval(point, &commitment, &witness), "expected verify_eval to fail for for point {:#?}, commitment {:#?}, and witness {:#?}, but it didn't", point, commitment, witness);
    }

    #[test]
    fn test_eval_basic() {
        const N: usize = 16;
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<N>(&mut rng);
        let lagrange_basis = compute_lagrange_basis(&params);

        let (prover, verifier) = test_participants(&params, &lagrange_basis.0, &lagrange_basis.1);

        let evals = random_evals::<N>(&mut rng);
        let commitment = prover.commit(&evals);

        let witness = prover.create_witness(&evals, 3);
        assert_verify_eval(&verifier, (3, evals.coeffs[3]), &commitment, &witness);

        let y_prime = random_field_elem_neq(&mut rng, evals.coeffs[3]);
        assert_verify_eval_fails(&verifier, (3, y_prime), &commitment, &witness);
    }

    #[test]
    fn test_eval_all() {
        const N: usize = 16;
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup::<N>(&mut rng);
        let lagrange_basis = compute_lagrange_basis(&params);

        let (mut prover, verifier) = test_participants(&params, &lagrange_basis.0, &lagrange_basis.1);

        let evals = random_evals::<N>(&mut rng);
        let commitment = prover.commit(&evals);
        let witness = prover.create_witness_all();
        assert!(verifier.verify_eval_all(&evals.coeffs, &commitment, &witness))
    }
}
