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

fn bench_verify_eval<const NUM_COEFFS: usize>(c: &mut Criterion) {
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

    let x: Scalar = rng.gen::<u64>().into();
    let y = polynomial.eval(x);
    let witness = prover.create_witness(&polynomial, (x, y)).unwrap();

    c.bench_function(
        format!("bench_verify_eval_coeff_form, degree {}", NUM_COEFFS - 1).as_str(),
        |b| {
            b.iter(|| {
                verifier.verify_eval(
                    black_box((x, y)),
                    black_box(&commitment),
                    black_box(&witness),
                )
            })
        },
    );
}

mod perf;

criterion_group!(
    name = verify_eval;
    config = Criterion::default().with_profiler(perf::FlamegraphProfiler::new(100));
    targets = bench_verify_eval<16>, bench_verify_eval<64>, bench_verify_eval<128>, bench_verify_eval<256>
);
criterion_main!(verify_eval);