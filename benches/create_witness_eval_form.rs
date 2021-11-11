use blstrs::Scalar;
use kzg::ft::EvaluationDomain;
use kzg::{eval_form::KZGProverEvalForm, setup, eval_form::compute_lagrange_basis, KZGParams};
use pairing::group::ff::Field;
use rand::{rngs::SmallRng, Rng, SeedableRng};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn test_setup(rng: &mut SmallRng, d: usize) -> KZGParams {
    let s: Scalar = rng.gen::<u64>().into();
    setup(s, d)
}

fn random_evals(rng: &mut SmallRng, d: usize) -> EvaluationDomain {
    let mut coeffs = vec![Scalar::zero(); d];

    for i in 0..d {
        coeffs[i] = rng.gen::<u64>().into();
    }

    EvaluationDomain::from_coeffs(coeffs).unwrap()
}

fn bench_create_witness<const NUM_COEFFS: usize>(c: &mut Criterion) {
    let mut rng = SmallRng::from_seed([42; 32]);
    let params = test_setup(&mut rng, NUM_COEFFS);
    let lagrange_basis = compute_lagrange_basis(&params);

    let evals = random_evals(&mut rng, NUM_COEFFS);
    let prover = KZGProverEvalForm::new(&params, lagrange_basis.0.as_slice(), lagrange_basis.1.as_slice());
    let _commitment = prover.commit(&evals);

    c.bench_function(
        format!("bench_create_witness_eval_form, degree {}", NUM_COEFFS - 1).as_str(),
        |b| b.iter(|| black_box(&prover).create_witness(black_box(&evals), black_box(5))),
    );

    c.bench_function(
        format!(
            "bench_create_witness_batched_eval_form, degree {}",
            NUM_COEFFS - 1
        )
        .as_str(),
        |b| {
            b.iter(|| {
                black_box(&prover)
                    .create_witness_all()
            })
        },
    );
}

criterion_group!(
    create_witness,
    bench_create_witness<16>,
    bench_create_witness<64>,
    bench_create_witness<128>,
    bench_create_witness<256>
);
criterion_main!(create_witness);
