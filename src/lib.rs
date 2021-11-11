use blstrs::{G1Affine, G1Projective, G2Projective, Scalar};
use pairing::group::{Curve, Group, prime::PrimeCurveAffine};
use thiserror::Error;

pub mod coeff_form;
pub mod eval_form;
pub mod ft;
pub mod polynomial;
pub mod utils;

/// parameters from tested setup
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde_support", derive(Serialize, Deserialize))]
pub struct KZGParams {
    /// g, g^alpha^1, g^alpha^2, ...
    gs: Vec<G1Projective>,
    /// h, h^alpha^1, h^alpha^2, ...
    hs: Vec<G2Projective>,
}

/// the commitment - "C" in the paper. It's a single group element
pub type KZGCommitment = G1Affine;
/// A witness for a single element - "w_i" in the paper. It's a group element.
pub type KZGWitness = G1Affine;

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

pub fn setup(s: Scalar, num_coeffs: usize) -> KZGParams {
    let mut gs = vec![G1Projective::generator(); num_coeffs];
    let mut hs = vec![G2Projective::generator(); num_coeffs];

    let mut curr = gs[0];
    for g in gs.iter_mut().skip(1) {
        *g = curr * s;
        curr = *g;
    }

    let mut curr = hs[0];
    for h in hs.iter_mut().skip(1) {
        *h = curr * s;
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
