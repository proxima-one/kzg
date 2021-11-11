use blstrs::Scalar;
use kzg::polynomial::Polynomial;
use kzg::{
    coeff_form::{KZGProver, KZGVerifier},
    setup, KZGParams,
};
use pairing::group::ff::Field;
use rand::{rngs::SmallRng, Rng, SeedableRng};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

pub fn csprng_setup<const MAX_COEFFS: usize>() -> KZGParams {
    let s: Scalar = rand::random::<u64>().into();
    setup(s, MAX_COEFFS)
}

fn bench_poly_verify<const NUM_COEFFS: usize>(c: &mut Criterion) {
    let params = csprng_setup::<NUM_COEFFS>();
    let mut rng = SmallRng::from_seed([42; 32]);
    let mut coeffs = vec![Scalar::zero(); NUM_COEFFS];
    for i in 0..NUM_COEFFS {
        coeffs[i] = rng.gen::<u64>().into();
    }
    let polynomial = Polynomial::new_from_coeffs(coeffs, NUM_COEFFS - 1);
    let prover = KZGProver::new(&params);
    let verifier = KZGVerifier::new(&params);
    let commitment = prover.commit(&polynomial);

    c.bench_function(
        format!("bench_verify_poly_coeff_form, degree {}", NUM_COEFFS - 1).as_str(),
        |b| {
            b.iter(|| {
                verifier.verify_poly(black_box(&commitment), black_box(&polynomial));
            })
        },
    );
}

criterion_group!(
    poly_verify,
    bench_poly_verify<16>,
    bench_poly_verify<64>,
    bench_poly_verify<128>,
    bench_poly_verify<256>
);
criterion_main!(poly_verify);
