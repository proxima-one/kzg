use bls12_381::Scalar;
use kzg::{polynomial::Polynomial, worker::Worker};
use pairing::group::ff::{Field, PrimeField};
use rand::{rngs::SmallRng, Rng, SeedableRng};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn random_polynomial<S: PrimeField>(rng: &mut SmallRng, n: usize) -> Polynomial<S> {
    let mut coeffs = vec![S::zero(); n];
    for i in 0..n {
        coeffs[i] = rng.gen::<u64>().into();
    }
    Polynomial::new(coeffs)
}

fn bench_poly_arithmetic<S: PrimeField, const NUM_COEFFS: usize>(c: &mut Criterion) {
    let mut rng = SmallRng::from_seed([NUM_COEFFS as u8; 32]);
    let f = random_polynomial::<S>(&mut rng, NUM_COEFFS);
    let g = random_polynomial::<S>(&mut rng, NUM_COEFFS);

    c.bench_function(
        format!("bench_add, degree {}", NUM_COEFFS - 1).as_str(),
        |b| {
            b.iter(|| black_box(f.clone()) + black_box(g.clone()));
        },
    );

    c.bench_function(
        format!("bench_mul_naive, degree {}", NUM_COEFFS - 1).as_str(),
        |b| {
            b.iter(|| black_box(f.clone()) * black_box(g.clone()));
        },
    );

    let worker = Worker::new();

    c.bench_function(
        format!("bench_mul_fft, degree {}", NUM_COEFFS - 1).as_str(),
        |b| {
            b.iter(|| black_box(f.clone()).fft_mul(black_box(g.clone()), &worker));
        },
    );

    let g = random_polynomial::<S>(&mut rng, 2);

    c.bench_function(
        format!("bench_long_division, degree {}", NUM_COEFFS - 1).as_str(),
        |b| b.iter(|| black_box(f.clone()).long_division(&black_box(g.clone()))),
    );

    let mut xs = vec![Scalar::zero(); NUM_COEFFS - 1];
    let mut ys = vec![Scalar::zero(); NUM_COEFFS - 1];
    for i in 0..xs.len() {
        xs[i] = Scalar::random(&mut rng);
        ys[i] = Scalar::random(&mut rng);
    }

    c.bench_function(
        format!("bench_interpolation, degree {}", NUM_COEFFS - 1).as_str(),
        |b| {
            b.iter(|| {
                Polynomial::lagrange_interpolation(
                    xs.as_slice(),
                    ys.as_slice(),
                )
            })
        },
    );
}

criterion_group!(poly_arithmetic, bench_poly_arithmetic<Scalar, 10>, bench_poly_arithmetic<Scalar, 50>, bench_poly_arithmetic<Scalar, 100>, bench_poly_arithmetic<Scalar, 200>);
criterion_main!(poly_arithmetic);
