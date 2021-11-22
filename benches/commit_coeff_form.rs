use blstrs::Scalar;
use kzg::polynomial::Polynomial;
use kzg::{coeff_form::KZGProver, setup, KZGParams};
use pairing::group::ff::Field;
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
    let prover = KZGProver::new(&params);

    c.bench_function(
        format!("bench_commit_coeff_form, degree {}", NUM_COEFFS - 1).as_str(),
        |b| {
            b.iter(|| {
                black_box(&prover).commit(black_box(&polynomial))
            })
        },
    );
}

mod perf;

criterion_group!(
    name = commit;
    config = Criterion::default().with_profiler(perf::FlamegraphProfiler::new(100));
    targets = bench_commit<16>, bench_commit<64>, bench_commit<128>, bench_commit<256>
);
criterion_main!(commit);
