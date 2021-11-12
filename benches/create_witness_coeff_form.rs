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

fn bench_create_witness<const NUM_COEFFS: usize>(c: &mut Criterion) {
    let params = csprng_setup::<NUM_COEFFS>();
    let mut rng = SmallRng::from_seed([42; 32]);
    let mut coeffs = vec![Scalar::zero(); NUM_COEFFS];
    for i in 0..NUM_COEFFS {
        coeffs[i] = rng.gen::<u64>().into();
    }
    let polynomial = Polynomial::new_from_coeffs(coeffs, NUM_COEFFS - 1);
    let prover = KZGProver::new(&params);
    let _commitment = prover.commit(&polynomial);

    let x: Scalar = Scalar::random(&mut rng);
    let y = polynomial.eval(x);
    
    c.bench_function(
        format!("bench_create_witness_coeff_form, degree {}", NUM_COEFFS - 1).as_str(),
        |b| b.iter(|| black_box(&prover).create_witness(black_box(&polynomial), black_box((x, y))).unwrap()),
    );

    let mut xs = Vec::with_capacity(NUM_COEFFS - 1);
    let mut ys = Vec::with_capacity(NUM_COEFFS - 1);
    for _ in 0..NUM_COEFFS - 1 {
        let x = Scalar::random(&mut rng);
        xs.push(x);
        ys.push(polynomial.eval(x));
    }
}

criterion_group!(
    create_witness,
    bench_create_witness<16>,
    bench_create_witness<64>,
    bench_create_witness<128>,
    bench_create_witness<256>
);
criterion_main!(create_witness);
