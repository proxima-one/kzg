use kzg::polynomial::Polynomial;
use kzg::{setup, KZGParams, KZGProver};
use pairing::{group::ff::Field};
use blstrs::Scalar;
use rand::{rngs::SmallRng, Rng, SeedableRng};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

pub fn csprng_setup<const MAX_COEFFS: usize>() -> KZGParams {
    let s: Scalar = rand::random::<u64>().into();
    setup(s, MAX_COEFFS)
}

fn bench_commit<const NUM_COEFFS: usize>(c: &mut Criterion) {
    let params = csprng_setup::<NUM_COEFFS>();
    let mut rng = SmallRng::from_seed([42; 32]);
    let mut coeffs = vec![Scalar::zero(); NUM_COEFFS];
    for i in 0..NUM_COEFFS {
        coeffs[i] = rng.gen::<u64>().into();
    }
    let polynomial = Polynomial::new(coeffs);

    c.bench_function(
        format!("bench_commit, degree {}", NUM_COEFFS - 1).as_str(),
        |b| {
            b.iter(|| {
                let mut prover = KZGProver::new(&params);
                prover.commit(black_box(polynomial.clone()))
            })
        },
    );
}

criterion_group!(commit, bench_commit<10>, bench_commit<50>, bench_commit<100>, bench_commit<200>);
criterion_main!(commit);
