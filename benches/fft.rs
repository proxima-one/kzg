use blstrs::Scalar;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use kzg::ft::EvaluationDomain;
use pairing::group::ff::Field;
use rand::{rngs::SmallRng, SeedableRng};

fn random_vector<const N: usize>(rng: &mut SmallRng) -> Vec<Scalar> {
    let mut coeffs = Vec::with_capacity(N);
    for _ in 0..N {
        coeffs.push(Scalar::random(&mut *rng))
    }
    coeffs
}

fn bench_fft<const N: usize>(c: &mut Criterion) {
    let mut rng = SmallRng::from_seed([N as u8; 32]);
    let coeffs = random_vector::<N>(&mut rng);
    let domain = EvaluationDomain::from_coeffs(coeffs).unwrap();

    c.bench_function(format!("bench_fft, size {}", N).as_str(), |b| {
        b.iter(|| black_box(domain.clone()).fft())
    });

    c.bench_function(format!("bench_ifft, size {}", N).as_str(), |b| {
        b.iter(|| black_box(domain.clone()).ifft())
    });
}


mod perf;
criterion_group!(
    name = fft;
    config = Criterion::default().with_profiler(perf::FlamegraphProfiler::new(100));
    targets = bench_fft<16>, bench_fft<64>, bench_fft<128>, bench_fft<256>, bench_fft<512>
);
criterion_main!(fft);
