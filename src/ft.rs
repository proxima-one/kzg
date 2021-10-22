use crate::polynomial::Polynomial;
use crate::worker::Worker;
use crate::KZGError;
use pairing::{
    group::ff::{Field, PrimeField},
    Engine,
};
use std::convert::TryFrom;
use std::ops::{AddAssign, MulAssign, SubAssign};

pub struct EvaluationDomain<S: PrimeField> {
    coeffs: Vec<S>,
    exp: u32,
    omega: S,
    omegainv: S,
    geninv: S,
    minv: S,
}

impl<S: PrimeField> From<EvaluationDomain<S>> for Polynomial<S> {
    fn from(domain: EvaluationDomain<S>) -> Polynomial<S> {
        Polynomial::new(domain.coeffs)
    }
}

impl<S: PrimeField> AsRef<[S]> for EvaluationDomain<S> {
    fn as_ref(&self) -> &[S] {
        &self.coeffs
    }
}

impl<S: PrimeField> AsMut<[S]> for EvaluationDomain<S> {
    fn as_mut(&mut self) -> &mut [S] {
        &mut self.coeffs
    }
}

impl<S: PrimeField> EvaluationDomain<S> {
    pub fn into_coeffs(self) -> Vec<S> {
        self.coeffs
    }

    pub fn from_coeffs(mut coeffs: Vec<S>) -> Result<EvaluationDomain<S>, KZGError> {
        // Compute the size of our evaluation domain
        let mut m = 1;
        let mut exp = 0;

        // TODO cache this in a lazy static
        while m < coeffs.len() {
            m *= 2;
            exp += 1;

            // The pairing-friendly curve may not be able to support
            // large enough (radix2) evaluation domains.
            if exp >= S::S {
                return Err(KZGError::PolynomialDegreeTooLarge);
            }
        }

        // Compute omega, the 2^exp primitive root of unity
        let mut omega = S::root_of_unity();
        for _ in exp..S::S {
            omega = omega.square();
        }

        // Extend the coeffs vector with zeroes if necessary
        coeffs.resize(m, S::zero());

        Ok(EvaluationDomain {
            coeffs,
            exp,
            omega,
            omegainv: omega.invert().unwrap(),
            geninv: S::multiplicative_generator().invert().unwrap(),
            minv: S::from(m as u64).invert().unwrap(),
        })
    }

    pub fn fft(&mut self, worker: &Worker) {
        best_fft(&mut self.coeffs, worker, &self.omega, self.exp);
    }

    pub fn ifft(&mut self, worker: &Worker) {
        best_fft(&mut self.coeffs, worker, &self.omegainv, self.exp);

        worker.scope(self.coeffs.len(), |scope, chunk| {
            let minv = self.minv;

            for v in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_scope| {
                    for v in v {
                        v.mul_assign(&minv);
                    }
                });
            }
        });
    }

    pub fn distribute_powers(&mut self, worker: &Worker, g: S) {
        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (i, v) in self.coeffs.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_scope| {
                    let mut u = g.pow_vartime(&[(i * chunk) as u64]);
                    for v in v.iter_mut() {
                        v.mul_assign(&u);
                        u.mul_assign(&g);
                    }
                });
            }
        });
    }

    pub fn coset_fft(&mut self, worker: &Worker) {
        self.distribute_powers(worker, S::multiplicative_generator());
        self.fft(worker);
    }

    pub fn icoset_fft(&mut self, worker: &Worker) {
        let geninv = self.geninv;

        self.ifft(worker);
        self.distribute_powers(worker, geninv);
    }

    /// This evaluates t(tau) for this domain, which is
    /// tau^m - 1 for these radix-2 domains.
    pub fn z(&self, tau: &S) -> S {
        let mut tmp = tau.pow_vartime(&[self.coeffs.len() as u64]);
        tmp.sub_assign(&S::one());

        tmp
    }

    /// The target polynomial is the zero polynomial in our
    /// evaluation domain, so we must perform division over
    /// a coset.
    pub fn divide_by_z_on_coset(&mut self, worker: &Worker) {
        let i = self.z(&S::multiplicative_generator()).invert().unwrap();

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for v in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_scope| {
                    for v in v {
                        v.mul_assign(&i);
                    }
                });
            }
        });
    }

    /// Perform O(n) multiplication of two polynomials in the domain.
    pub fn mul_assign(&mut self, worker: &Worker, other: &EvaluationDomain<S>) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_scope| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        a.mul_assign(b);
                    }
                });
            }
        });
    }

    /// Perform O(n) subtraction of one polynomial from another in the domain.
    pub fn sub_assign(&mut self, worker: &Worker, other: &EvaluationDomain<S>) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_scope| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        a.sub_assign(b);
                    }
                });
            }
        });
    }
}

fn best_fft<S: PrimeField>(a: &mut [S], worker: &Worker, omega: &S, log_n: u32) {
    let log_cpus = worker.log_num_threads();

    if log_n <= log_cpus {
        serial_fft(a, omega, log_n);
    } else {
        parallel_fft(a, worker, omega, log_n, log_cpus);
    }
}

#[allow(clippy::many_single_char_names)]
fn serial_fft<S: PrimeField>(a: &mut [S], omega: &S, log_n: u32) {
    fn bitreverse(mut n: u32, l: u32) -> u32 {
        let mut r = 0;
        for _ in 0..l {
            r = (r << 1) | (n & 1);
            n >>= 1;
        }
        r
    }

    let n = a.len() as u32;
    assert_eq!(n, 1 << log_n);

    for k in 0..n {
        let rk = bitreverse(k, log_n);
        if k < rk {
            a.swap(rk as usize, k as usize);
        }
    }

    let mut m = 1;
    for _ in 0..log_n {
        let w_m = omega.pow_vartime(&[u64::from(n / (2 * m))]);

        let mut k = 0;
        while k < n {
            let mut w = S::one();
            for j in 0..m {
                let mut t = a[(k + j + m) as usize];
                t.mul_assign(&w);
                let mut tmp = a[(k + j) as usize];
                tmp.sub_assign(&t);
                a[(k + j + m) as usize] = tmp;
                a[(k + j) as usize].add_assign(&t);
                w.mul_assign(&w_m);
            }

            k += 2 * m;
        }

        m *= 2;
    }
}

fn parallel_fft<S: PrimeField>(a: &mut [S], worker: &Worker, omega: &S, log_n: u32, log_cpus: u32) {
    assert!(log_n >= log_cpus);

    let num_cpus = 1 << log_cpus;
    let log_new_n = log_n - log_cpus;
    let mut tmp = vec![vec![S::zero(); 1 << log_new_n]; num_cpus];
    let new_omega = omega.pow_vartime(&[num_cpus as u64]);

    worker.scope(0, |scope, _| {
        let a = &*a;

        for (j, tmp) in tmp.iter_mut().enumerate() {
            scope.spawn(move |_scope| {
                // Shuffle into a sub-FFT
                let omega_j = omega.pow_vartime(&[j as u64]);
                let omega_step = omega.pow_vartime(&[(j as u64) << log_new_n]);

                let mut elt = S::one();
                for (i, tmp) in tmp.iter_mut().enumerate() {
                    for s in 0..num_cpus {
                        let idx = (i + (s << log_new_n)) % (1 << log_n);
                        let mut t = a[idx];
                        t.mul_assign(&elt);
                        tmp.add_assign(&t);
                        elt.mul_assign(&omega_step);
                    }
                    elt.mul_assign(&omega_j);
                }

                // Perform sub-FFT
                serial_fft(tmp, &new_omega, log_new_n);
            });
        }
    });

    // TODO: does this hurt or help?
    worker.scope(a.len(), |scope, chunk| {
        let tmp = &tmp;

        for (idx, a) in a.chunks_mut(chunk).enumerate() {
            scope.spawn(move |_scope| {
                let mut idx = idx * chunk;
                let mask = (1 << log_cpus) - 1;
                for a in a {
                    *a = tmp[idx & mask][idx >> log_cpus];
                    idx += 1;
                }
            });
        }
    });
}

// Test multiplying various (low degree) polynomials together and
// comparing with naive evaluations.
#[test]
fn polynomial_arith() {
    use bls12_381::Scalar as Fr;
    use rand::RngCore;

    fn test_mul<S: PrimeField, R: RngCore>(mut rng: &mut R) {
        let worker = Worker::new();

        for coeffs_a in vec![1, 5, 10, 50] {
            for coeffs_b in vec![1, 5, 10, 50] {
                let a: Vec<_> = (0..coeffs_a).map(|_| S::random(&mut rng)).collect();
                let b: Vec<_> = (0..coeffs_b).map(|_| S::random(&mut rng)).collect();

                let a = Polynomial::new_from_coeffs(a, coeffs_a - 1);
                let b = Polynomial::new_from_coeffs(b, coeffs_b - 1);

                // naive evaluation
                let naive = a.clone() * b.clone();
                let fft = a.fft_mul(&b, &worker);

                assert!(naive == fft);
            }
        }
    }

    let rng = &mut rand::thread_rng();

    test_mul::<Fr, _>(rng);
}

#[test]
fn fft_composition() {
    use bls12_381::Scalar as Fr;
    use rand::RngCore;

    fn test_comp<S: PrimeField, R: RngCore>(mut rng: &mut R) {
        let worker = Worker::new();

        for coeffs in 0..10 {
            let coeffs = 1 << coeffs;

            let mut v = vec![];
            for _ in 0..coeffs {
                v.push(S::random(&mut rng));
            }

            let mut domain = EvaluationDomain::from_coeffs(v.clone()).unwrap();
            domain.ifft(&worker);
            domain.fft(&worker);
            assert!(v == domain.coeffs);
            domain.fft(&worker);
            domain.ifft(&worker);
            assert!(v == domain.coeffs);
            domain.icoset_fft(&worker);
            domain.coset_fft(&worker);
            assert!(v == domain.coeffs);
            domain.coset_fft(&worker);
            domain.icoset_fft(&worker);
            assert!(v == domain.coeffs);
        }
    }

    let rng = &mut rand::thread_rng();

    test_comp::<Fr, _>(rng);
}

#[test]
fn parallel_fft_consistency() {
    use bls12_381::Scalar as Fr;
    use rand::RngCore;
    use std::cmp::min;

    fn test_consistency<S: PrimeField, R: RngCore>(mut rng: &mut R) {
        let worker = Worker::new();

        for _ in 0..5 {
            for log_d in 0..10 {
                let d = 1 << log_d;

                let v1 = (0..d).map(|_| S::random(&mut rng)).collect::<Vec<_>>();
                let mut v1 = EvaluationDomain::from_coeffs(v1).unwrap();
                let mut v2 = EvaluationDomain::from_coeffs(v1.coeffs.clone()).unwrap();

                for log_cpus in log_d..min(log_d + 1, 3) {
                    parallel_fft(&mut v1.coeffs, &worker, &v1.omega, log_d, log_cpus);
                    serial_fft(&mut v2.coeffs, &v2.omega, log_d);

                    assert!(v1.coeffs == v2.coeffs);
                }
            }
        }
    }

    let rng = &mut rand::thread_rng();

    test_consistency::<Fr, _>(rng);
}
