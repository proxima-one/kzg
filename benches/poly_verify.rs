use bls12_381::Bls12;
use kzg::polynomial::Polynomial;
use kzg::{setup, KZGParams, KZGProver, KZGVerifier};
use pairing::{group::ff::Field, Engine};
use rand::{rngs::SmallRng, Rng, SeedableRng};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

pub fn csprng_setup<E: Engine, const MAX_COEFFS: usize>() -> KZGParams<E> {
    let s: E::Fr = rand::random::<u64>().into();
    setup(s, MAX_COEFFS)
}

fn bench_poly_verify<E: Engine, const NUM_COEFFS: usize>(c: &mut Criterion) {
    let params = csprng_setup::<E, NUM_COEFFS>();
    let mut rng = SmallRng::from_seed([42; 32]);
    let mut coeffs = vec![E::Fr::zero(); NUM_COEFFS];
    for i in 0..NUM_COEFFS {
        coeffs[i] = rng.gen::<u64>().into();
    }
    let polynomial = Polynomial::new_from_coeffs(coeffs, NUM_COEFFS - 1);
    let mut prover = KZGProver::new(&params);
    let verifier = KZGVerifier::new(&params);
    let commitment = prover.commit(polynomial.clone());

    c.bench_function(
        format!("bench_verify_poly, degree {}", NUM_COEFFS - 1).as_str(),
        |b| {
            b.iter(|| {
                verifier.verify_poly(black_box(&commitment), black_box(&polynomial));
            })
        },
    );
}

criterion_group!(poly_verify, bench_poly_verify<Bls12, 10>, bench_poly_verify<Bls12, 50>, bench_poly_verify<Bls12, 100>, bench_poly_verify<Bls12, 200>);
criterion_main!(poly_verify);
