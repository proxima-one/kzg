use blstrs::{pairing, G1Affine, G1Projective, G2Affine, G2Projective, Scalar};
use pairing::group::ff::PrimeField;
use pairing::group::{ff::Field, Group, prime::PrimeCurveAffine, Curve};
use std::fmt::Debug;

#[cfg(feature = "serde_support")]
use serde::{Deserialize, Serialize};

use crate::ft::EvaluationDomain;
use crate::polynomial::Polynomial;
use crate::{KZGCommitment, KZGError, KZGParams, KZGWitness};

// A witness for a several elements - "w_B" in the paper. It's a single group element plus a polynomial
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde_support", derive(Serialize, Deserialize))]
pub struct KZGBatchWitnessEvalForm {
    r: EvaluationDomain,
    w: G1Affine,
}

impl KZGBatchWitnessEvalForm {
    pub fn elem(&self) -> G1Affine {
        self.w
    }

    pub fn elem_ref(&self) -> &G1Affine {
        &self.w
    }

    pub fn polynomial(&self) -> &EvaluationDomain {
        &self.r
    }

    pub fn new(r: EvaluationDomain, w: G1Affine) -> Self {
        KZGBatchWitnessEvalForm { r, w }
    }
}

#[derive(Debug, Clone)]
pub struct KZGProverEvalForm<'params> {
    parameters: &'params KZGParams,
    lagrange_basis_g: &'params [G1Projective],
    lagrange_basis_h: &'params [G2Projective],
    d: usize,
    exp: u32,
    omega: Scalar,
}

#[derive(Debug, Clone)]
pub struct KZGVerifierEvalForm<'params> {
    d: usize,
    exp: u32,
    omega: Scalar,
    parameters: &'params KZGParams,
    lagrange_basis_g: &'params [G1Projective],
    lagrange_basis_h: &'params [G2Projective],
}

fn div_by_omega_i(evals: &EvaluationDomain, m: usize) -> EvaluationDomain {
    let mut coeffs = Vec::with_capacity(evals.d);
    let omega_m = evals.omega.pow_vartime(&[m as u64]);
    for (j, &f) in evals.coeffs.iter().enumerate() {
        if j == m {
            let mut qm = Scalar::zero();
            let am = Scalar::from(evals.d as u64) * omega_m.invert().unwrap();

            for i in (0..evals.d).filter(|&i| i != m) {
                let omega_i = evals.omega.pow_vartime(&[i as u64]);
                let ai = Scalar::from(evals.d as u64) * omega_i.invert().unwrap();

                let mut term = evals.coeffs[i];
                term *= am * ai.invert().unwrap();
                term *= (omega_m - omega_i).invert().unwrap();
                qm += term;
            }

            coeffs.push(qm);
        } else {
            let omega_j = evals.omega.pow_vartime(&[j as u64]);
            coeffs.push(f * (omega_j - omega_m).invert().unwrap())
        }
    }

    evals.clone_with_different_coeffs(coeffs)
}

impl<'params> KZGProverEvalForm<'params> {
    /// initializes `polynomial` to zero polynomial
    pub fn new(
        parameters: &'params KZGParams,
        lagrange_basis_g: &'params [G1Projective],
        lagrange_basis_h: &'params [G2Projective],
    ) -> Self {
        let (d, exp, omega) = EvaluationDomain::compute_omega(parameters.gs.len()).unwrap();
        Self {
            parameters,
            lagrange_basis_g,
            lagrange_basis_h,
            d,
            exp,
            omega,
        }
    }

    pub fn parameters(&self) -> &'params KZGParams {
        self.parameters
    }

    pub fn degree(&self) -> usize {
        self.d
    }

    pub fn omega(&self) -> Scalar {
        self.omega
    }

    pub fn commit(&self, evals: &EvaluationDomain) -> KZGCommitment {
        assert!(self.d == evals.d);

        let gs = &self.lagrange_basis_g[..evals.len()];
        let commitment = G1Projective::multi_exp(gs, evals.as_ref());

        let commitment = commitment.to_affine();
        commitment
    }

    pub fn create_witness(&self, evals: &EvaluationDomain, i: usize) -> KZGWitness {
        let y = evals.coeffs[i];
        let numerator = EvaluationDomain {
            coeffs: evals.coeffs.iter().map(|c| c - y).collect(),
            ..*evals
        };

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

impl<'params> KZGVerifierEvalForm<'params> {
    pub fn new(parameters: &'params KZGParams, lagrange_basis_g: &'params [G1Projective], lagrange_basis_h: &'params [G2Projective]) -> Self {
        let (d, exp, omega) = EvaluationDomain::compute_omega(parameters.gs.len()).unwrap();
        KZGVerifierEvalForm {
            parameters,
            d,
            exp,
            omega,
            lagrange_basis_g,
            lagrange_basis_h
        }
    }

    pub fn verify_poly(&self, commitment: &KZGCommitment, evals: &EvaluationDomain) -> bool {
        let mut evals = evals.clone();
        evals.ifft();
        let polynomial: Polynomial = evals.into();

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
        ys: &[Scalar],
        commitment: &KZGCommitment,
        witness: &KZGWitness,
    ) -> bool {
        let mut z = EvaluationDomain::new(vec![Scalar::zero(); self.d], self.d, self.exp, self.omega);
        z.coeffs[0] = -Scalar::one();
        z.coeffs[self.d - 1] = Scalar::one();

        let hs = &self.lagrange_basis_h[..z.len()];
        let hz = G2Projective::multi_exp(hs, z.coeffs.as_slice());

        let r = EvaluationDomain::new(ys.to_vec(), self.d, self.exp, self.omega);

        let gs = &self.lagrange_basis_g[..r.len()];
        let gr = G1Projective::multi_exp(gs, r.coeffs.as_slice());

        let lhs = pairing(witness, &hz.to_affine());
        let rhs = pairing(
            &(commitment.to_curve() - gr).to_affine(),
            &self.parameters.hs[0].to_affine(),
        );

        lhs == rhs
    }
}

pub fn compute_lagrange_basis_and_polynomials(params: &KZGParams) -> (Vec<G1Projective>, Vec<G2Projective>, Vec<Polynomial>) {
    let d = params.gs.len();
    assert!(d & (d - 1) == 0);

    let (d, _, omega) = EvaluationDomain::compute_omega(params.gs.len()).unwrap();
    let mut gs = Vec::with_capacity(d);
    let mut hs = Vec::with_capacity(d);
    let mut ls = Vec::with_capacity(d);

    for i in 0..d {
        let xi = omega.pow_vartime(&[i as u64]);
        let mut l = Polynomial::new_monic_of_degree(0);
        for j in (0..d).filter(|&j| j != i) {
            let xj = omega.pow_vartime(&[j as u64]);
            l = l.best_mul(&Polynomial::new(vec![-xj, Scalar::one()]));
            l = l.scalar_multiplication((xi - xj).invert().unwrap());
        }


        let coeffs = l.coeffs.clone();
        let g = &params.gs[..coeffs.len()];
        gs.push(G1Projective::multi_exp(g, coeffs.as_slice()));

        let h = &params.hs[..coeffs.len()];
        hs.push(G2Projective::multi_exp(h, coeffs.as_slice()));

        ls.push(l);
    }

    (gs, hs, ls)
}

/// params.gs.len() must be a power of two.
pub fn compute_lagrange_basis(params: &KZGParams) -> (Vec<G1Projective>, Vec<G2Projective>) {
    let d = params.gs.len();
    assert!(d & (d - 1) == 0);

    let (d, _, omega) = EvaluationDomain::compute_omega(params.gs.len()).unwrap();
    let mut gs = Vec::with_capacity(d);
    let mut hs = Vec::with_capacity(d);

    for i in 0..d {
        let xi = omega.pow_vartime(&[i as u64]);
        let mut l = Polynomial::new_monic_of_degree(0);
        for j in (0..d).filter(|&j| j != i) {
            let xj = omega.pow_vartime(&[j as u64]);
            l = l.best_mul(&Polynomial::new(vec![-xj, Scalar::one()]));
            l = l.scalar_multiplication((xi - xj).invert().unwrap());
        }

        let coeffs = l.coeffs();
        let g = &params.gs[..coeffs.len()];
        gs.push(G1Projective::multi_exp(g, coeffs.as_slice()));

        let h = &params.hs[..coeffs.len()];
        hs.push(G2Projective::multi_exp(h, coeffs.as_slice()));
    }

    (gs, hs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::setup;
    use rand::{rngs::SmallRng, Rng, SeedableRng};

    const RNG_SEED: [u8; 32] = [69; 32];

    fn test_setup(rng: &mut SmallRng, d: usize) -> KZGParams {
        assert!(d & (d - 1) == 0);

        let s: Scalar = rng.gen::<u64>().into();
        setup(s, d)
    }

    fn test_participants<'params>(
        params: &'params KZGParams,
        lagrange_basis_g: &'params [G1Projective],
        lagrange_basis_h: &'params [G2Projective],
    ) -> (KZGProverEvalForm<'params>, KZGVerifierEvalForm<'params>) {
        let prover = KZGProverEvalForm::new(params, lagrange_basis_g, lagrange_basis_h);
        let verifier = KZGVerifierEvalForm::new(params, lagrange_basis_g, lagrange_basis_h);

        (prover, verifier)
    }

    fn random_evals(rng: &mut SmallRng, d: usize) -> EvaluationDomain {
        let mut coeffs = vec![Scalar::zero(); d];

        for i in 0..d {
            coeffs[i] = rng.gen::<u64>().into();
        }

        EvaluationDomain::from_coeffs(coeffs).unwrap()
    }

    #[test]
    fn test_div_by_omega_i() {
        const N: usize = 10;
        let (d, exp, omega) = EvaluationDomain::compute_omega(N).unwrap();
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let omega_3 = omega.pow_vartime(&[3 as u64]);

        let mut top = Polynomial::new(random_evals(&mut rng, d).coeffs);
        let y = top.eval(omega_3);
        top.coeffs[0] -= y;
        let bottom = Polynomial::new(vec![-omega_3, Scalar::one()]);
        let (naive, rem) = top.long_division(&bottom);
        assert!(rem.is_none());

        let mut top = EvaluationDomain::new(top.coeffs, d, exp, omega);
        top.fft();
        let mut smart = div_by_omega_i(&top, 3);
        smart.ifft();
        let smart: Polynomial = smart.into();

        assert_eq!(smart, naive);
    }

    #[test]
    fn test_subtract_by_scalar() {
        const N: usize = 10;
        let (d, _, _) = EvaluationDomain::compute_omega(N).unwrap();
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let mut evals = random_evals(&mut rng, d);
        let subend = evals.clone_with_different_coeffs(vec![2.into(); d]);

        let mut naive = evals.clone();
        naive.ifft();
        let mut naive: Polynomial = naive.into();
        naive.coeffs[0] -= Scalar::from(2);

        evals.sub_assign(&subend);
        evals.ifft();
        let smart: Polynomial = evals.into();
        assert_eq!(smart, naive);
    }

    fn assert_verify_poly(
        verifier: &KZGVerifierEvalForm,
        commitment: &KZGCommitment,
        evals: &EvaluationDomain,
    ) {
        assert!(
            verifier.verify_poly(&commitment, &evals),
            "verify_poly failed for commitment {:#?} and polynomial {:#?}",
            commitment,
            evals
        );
    }

    fn assert_verify_poly_fails(
        verifier: &KZGVerifierEvalForm,
        commitment: &KZGCommitment,
        evals: &EvaluationDomain,
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
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup(&mut rng, 16);
        let lagrange_basis = compute_lagrange_basis(&params);

        let (mut prover, verifier) = test_participants(&params, lagrange_basis.0.as_slice(), lagrange_basis.1.as_slice());

        let evals = random_evals(&mut rng, prover.d);
        let commitment = prover.commit(&evals);

        assert_verify_poly(&verifier, &commitment, &evals);
        assert_verify_poly_fails(&verifier, &commitment, &random_evals(&mut rng, prover.d));
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
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup(&mut rng, 8);
        let lagrange_basis = compute_lagrange_basis(&params);

        let (prover, verifier) = test_participants(&params, lagrange_basis.0.as_slice(), lagrange_basis.1.as_slice());

        let evals = random_evals(&mut rng, prover.d);
        let commitment = prover.commit(&evals);

        let mut modified_evals = evals.clone();
        let new_coeff = random_field_elem_neq(&mut rng, modified_evals.coeffs[2]);
        modified_evals.coeffs[2] = new_coeff;

        assert_verify_poly(&verifier, &commitment, &evals);
        assert_verify_poly_fails(&verifier, &commitment, &modified_evals);
    }

    fn assert_verify_eval(
        verifier: &KZGVerifierEvalForm,
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

    fn assert_verify_eval_fails(
        verifier: &KZGVerifierEvalForm,
        point: (usize, Scalar),
        commitment: &KZGCommitment,
        witness: &KZGWitness,
    ) {
        assert!(!verifier.verify_eval(point, &commitment, &witness), "expected verify_eval to fail for for point {:#?}, commitment {:#?}, and witness {:#?}, but it didn't", point, commitment, witness);
    }

    #[test]
    fn test_eval_basic() {
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup(&mut rng, 16);
        let lagrange_basis = compute_lagrange_basis(&params);

        let (prover, verifier) = test_participants(&params, lagrange_basis.0.as_slice(), lagrange_basis.1.as_slice());

        let evals = random_evals(&mut rng, prover.d);
        let commitment = prover.commit(&evals);

        let witness = prover.create_witness(&evals, 3);
        assert_verify_eval(&verifier, (3, evals.coeffs[3]), &commitment, &witness);

        let y_prime = random_field_elem_neq(&mut rng, evals.coeffs[3]);
        assert_verify_eval_fails(&verifier, (3, y_prime), &commitment, &witness);
    }

    #[test]
    fn test_eval_all() {
        let mut rng = SmallRng::from_seed(RNG_SEED);
        let params = test_setup(&mut rng, 16);
        let lagrange_basis = compute_lagrange_basis(&params);

        let (mut prover, verifier) = test_participants(&params, lagrange_basis.0.as_slice(), lagrange_basis.1.as_slice());

        let evals = random_evals(&mut rng, prover.d);
        let commitment = prover.commit(&evals);
        let witness = prover.create_witness_all();
        assert!(verifier.verify_eval_all(evals.coeffs.as_ref(), &commitment, &witness))
    }
}
