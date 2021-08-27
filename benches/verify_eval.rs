use bls12_381::Bls12;
use kzg::polynomial::Polynomial;
use kzg::{setup, KZGParams, KZGProver, KZGVerifier};
use pairing::{group::ff::Field, Engine};
use rand::{rngs::SmallRng, Rng, SeedableRng};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

pub fn csprng_setup<E: Engine, const MAX_COEFFS: usize>() -> KZGParams<E, MAX_COEFFS> {
    let s: E::Fr = rand::random::<u64>().into();
    setup(s)
}

fn bench_verify_eval<E: Engine, const NUM_COEFFS: usize>(c: &mut Criterion) {
    let params = csprng_setup::<E, NUM_COEFFS>();
    let mut rng = SmallRng::from_seed([42; 32]);
    let mut coeffs = [E::Fr::zero(); NUM_COEFFS];
    for i in 0..NUM_COEFFS {
        coeffs[i] = rng.gen::<u64>().into();
    }
    let polynomial: Polynomial<E, NUM_COEFFS> = Polynomial::new_from_coeffs(coeffs, NUM_COEFFS - 1);
    let mut prover = KZGProver::new(&params);
	let verifier = KZGVerifier::new(&params);
	let commitment = prover.commit(polynomial.clone());

	let x: E::Fr = rng.gen::<u64>().into();
	let y = polynomial.eval(x);
	let witness = prover.create_witness((x, y)).unwrap();

    c.bench_function(
        format!("bench_verify_eval, degree {}", NUM_COEFFS - 1).as_str(),
        |b| {
            b.iter(|| {
               verifier.verify_eval(black_box((x, y)), black_box(&commitment), black_box(&witness))
            })
        },
    );
}

criterion_group!(verify_eval, bench_verify_eval<Bls12, 10>, bench_verify_eval<Bls12, 50>, bench_verify_eval<Bls12, 100>, bench_verify_eval<Bls12, 200>, bench_verify_eval<Bls12, 500>);
criterion_main!(verify_eval);
