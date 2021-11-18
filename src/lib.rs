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
pub struct KZGParams<const N: usize> {
    /// g, g^alpha^1, g^alpha^2, ...
    gs: [G1Projective; N],
    /// h, h^alpha^1, h^alpha^2, ...
    hs: [G2Projective; N]
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

pub fn setup<const N: usize>(s: Scalar) -> KZGParams<N> {
    assert!(N & (N - 1) == 0);
    let mut gs = [G1Projective::generator(); N];
    let mut hs = [G2Projective::generator(); N];

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
pub fn csprng_setup<const N: usize>() -> KZGParams<N> {
    let s: Scalar = random::<u64>().into();
    setup::<N>(s)
}
