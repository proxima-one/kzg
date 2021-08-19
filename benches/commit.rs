use bls12_381::Bls12;
use kzg::polynomial::Polynomial;
use kzg::{setup, KZGParams, KZGProver};
use pairing::{group::ff::Field, Engine};
use rand::{rngs::SmallRng, Rng, SeedableRng};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

pub fn csprng_setup<E: Engine, const DEGREE_LIMIT: usize>() -> KZGParams<E, DEGREE_LIMIT> {
    let s: E::Fr = rand::random::<u64>().into();
    setup(s)
}

fn bench_commit<E: Engine, const NUM_COEFFS: usize>(c: &mut Criterion) {
    let params = csprng_setup::<E, NUM_COEFFS>();
    let mut rng = SmallRng::from_seed([42; 32]);
    let mut coeffs = [E::Fr::zero(); NUM_COEFFS];
    for i in 0..NUM_COEFFS {
        coeffs[i] = rng.gen::<u64>().into();
    }
    let polynomial = Polynomial::new_from_coeffs(coeffs, NUM_COEFFS - 1);

    c.bench_function(
        format!("bench_commit, degree {}", NUM_COEFFS - 1).as_str(),
        |b| {
            b.iter(|| {
                let mut prover = KZGProver::new(params.clone());
                prover.commit(black_box(black_box(polynomial.clone())))
            })
        },
    );
}

criterion_group!(commit, bench_commit<Bls12, 10>, bench_commit<Bls12, 50>, bench_commit<Bls12, 100>, bench_commit<Bls12, 200>, bench_commit<Bls12, 500>);
criterion_main!(commit);
