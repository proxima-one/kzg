## kzg

**This code has not beed audited - use it at your own risk.**

`kzg` is a simple implementation of the [Kate-Zaverucha-Goldberg polynomial commitment scheme](https://www.iacr.org/archive/asiacrypt2010/6477178/6477178.pdf) over the [`zkcrypto`](https://github.com/zkcrypto) ecosystem's primitives, mainly their [`pairing`](https://github.com/zkcrypto/pairing) abstraction.

`kzg` implements the "simple" variant described in the paper as "DL", including batched openings.

### Author's Note

I wrote this mostly to learn and partly because the [`arkworks-polycommit`](https://github.com/arkworks-rs/poly-commit/tree/master/src) is hard to use and only implements the pederson variant of KZG, which is unnecessary for many use cases, in particular vector commitment schemes that don't care about the unconditional hiding property the pederson variant of KZG provides like [this](https://ethresear.ch/t/open-problem-ideal-vector-commitment/7421).
