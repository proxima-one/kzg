use blstrs::Scalar;
use kzg::polynomial::Polynomial;
use pairing::group::ff::Field;
use rand::{rngs::SmallRng, Rng, SeedableRng};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn random_polynomial(rng: &mut SmallRng, n: usize) -> Polynomial {
    let mut coeffs = vec![Scalar::zero(); n];
    for i in 0..n {
        coeffs[i] = rng.gen::<u64>().into();
    }
    Polynomial::new(coeffs)
}

fn bench_poly_arithmetic<const NUM_COEFFS: usize>(c: &mut Criterion) {
    let mut rng = SmallRng::from_seed([NUM_COEFFS as u8; 32]);
    let f = random_polynomial(&mut rng, NUM_COEFFS);
    let g = random_polynomial(&mut rng, NUM_COEFFS);

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

    c.bench_function(
        format!("bench_mul_fft, degree {}", NUM_COEFFS - 1).as_str(),
        |b| {
            b.iter(|| black_box(f.clone()).fft_mul(black_box(&g)));
        },
    );

    let g = random_polynomial(&mut rng, 2);

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
        |b| b.iter(|| Polynomial::lagrange_interpolation(xs.as_slice(), ys.as_slice())),
    );
}

criterion_group!(
    poly_arithmetic,
    bench_poly_arithmetic<16>,
    bench_poly_arithmetic<64>,
    bench_poly_arithmetic<128>,
    bench_poly_arithmetic<256>,
    bench_poly_arithmetic<512>
);
criterion_main!(poly_arithmetic);
