use pairing::group::ff::PrimeField;

// fast 64-bit log
// copypasta from https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
const LOG_TABLE: [u64; 64] = [
    63, 0, 58, 1, 59, 47, 53, 2, 60, 39, 48, 27, 54, 33, 42, 3, 61, 51, 37, 40, 49, 18, 28, 20, 55,
    30, 34, 11, 43, 14, 22, 4, 62, 57, 46, 52, 38, 26, 32, 41, 50, 36, 17, 19, 29, 10, 13, 21, 56,
    45, 25, 31, 35, 16, 9, 12, 44, 24, 15, 8, 23, 7, 6, 5,
];

pub fn log2(mut x: u64) -> u64 {
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    LOG_TABLE[(((x - (x >> 1)) * 0x07EDD5E59A4E28C2) >> 58) as usize]
}

pub fn log2_ceil(x: u64) -> u64 {
    let n = log2(x);
    // if x is a power of two, n is ceiling, otherwise it's n + 1
    if x & (x - 1) == 0 {
        n
    } else {
        n + 1
    }
}

pub fn pad_to_power_of_two<S: PrimeField>(xs: &[S]) -> Vec<S> {
    let n = 1 << log2_ceil(xs.len() as u64) as usize;
    let mut xs: Vec<S> = xs.to_vec();
    if xs.len() != n {
        xs.resize(n, S::zero())
    }
    xs
}

#[cfg(feature = "parallel")]
pub fn chunk_by_num_threads(size: usize) -> usize {
    let num_threads = rayon::current_num_threads();
    if size < num_threads {
        1
    } else {
        size / num_threads
    }
}
