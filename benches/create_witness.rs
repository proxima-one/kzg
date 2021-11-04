use bls12_381::Bls12;
use kzg::polynomial::Polynomial;
use kzg::{setup, KZGParams, KZGProver};
use pairing::{group::ff::Field, Engine};
use rand::{rngs::SmallRng, Rng, SeedableRng};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

pub fn csprng_setup<E: Engine, const MAX_COEFFS: usize>() -> KZGParams<E> {
    let s: E::Fr = rand::random::<u64>().into();
    setup(s, MAX_COEFFS)
}

fn bench_create_witness<E: Engine, const NUM_COEFFS: usize>(c: &mut Criterion) {
    let params = csprng_setup::<E, NUM_COEFFS>();
    let mut rng = SmallRng::from_seed([42; 32]);
    let mut coeffs = vec![E::Fr::zero(); NUM_COEFFS];
    for i in 0..NUM_COEFFS {
        coeffs[i] = rng.gen::<u64>().into();
    }
    let polynomial = Polynomial::new_from_coeffs(coeffs, NUM_COEFFS - 1);
    let mut prover = KZGProver::new(&params);
    let _commitment = prover.commit(polynomial.clone());

    let x: E::Fr = E::Fr::random(&mut rng);
    let y = polynomial.eval(x);

    c.bench_function(
        format!("bench_create_witness, degree {}", NUM_COEFFS - 1).as_str(),
        |b| b.iter(|| prover.create_witness(black_box((x, y))).unwrap()),
    );

    let mut xs = Vec::with_capacity(NUM_COEFFS - 1);
    let mut ys = Vec::with_capacity(NUM_COEFFS - 1);
    for _ in 0..NUM_COEFFS - 1 {
        let x = E::Fr::random(&mut rng);
        xs.push(x);
        ys.push(polynomial.eval(x));
    }

    c.bench_function(
        format!("bench_create_witness_batched, degree {}", NUM_COEFFS - 1).as_str(),
        |b| b.iter(|| prover.create_witness_batched(black_box(xs.as_slice()), black_box(ys.as_slice())).unwrap()),
    );
}

criterion_group!(create_witness, bench_create_witness<Bls12, 10>, bench_create_witness<Bls12, 50>, bench_create_witness<Bls12, 100>, bench_create_witness<Bls12, 200>);
criterion_main!(create_witness);
