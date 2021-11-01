use pairing::{
    group::{GroupEncoding, ff::PrimeField},
    Engine,
};
use std::{fmt, ops::Deref, marker::PhantomData};

#[cfg(feature = "serde_support")]
use serde::{Serialize, Deserialize, Serializer, Deserializer, de::Visitor, de::Error as SerdeError, de::Unexpected};

// TODO: why is this PartialEq derive not working?
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd)]
pub struct G1Affine<E: Engine>(E::G1Affine);

impl<E: Engine> Copy for G1Affine<E> {}

impl<E: Engine> AsRef<E::G1Affine> for G1Affine<E> {
    fn as_ref(&self) -> &E::G1Affine {
        &self.0
    }
}

impl<E: Engine> AsMut<E::G1Affine> for G1Affine<E> {
    fn as_mut(&mut self) -> &mut E::G1Affine {
        &mut self.0
    }
}

impl<E: Engine> Deref for G1Affine<E> {
    type Target = E::G1Affine;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<E: Engine> G1Affine<E> {
    pub fn into_inner(self) -> E::G1Affine {
        self.0
    }

    pub fn inner(&self) -> &E::G1Affine {
        &self.0
    }

    pub fn from_inner(inner: E::G1Affine) -> Self {
        Self(inner)
    }
}

#[cfg(feature = "serde_support")]
impl<E: Engine> Serialize for G1Affine<E> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer
    {
        serializer.serialize_bytes(self.0.to_bytes().as_ref())
    }
}

#[cfg(feature = "serde_support")]
struct G1AffineVisitor<E>(PhantomData<E>);

#[cfg(feature = "serde_support")]
impl<'de, E: Engine> Visitor<'de> for G1AffineVisitor<E> {
    type Value = E::G1Affine;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        write!(formatter, "canonical bytes representation of an affine G1 Element")
    }

    fn visit_bytes<R>(self, v: &[u8]) -> Result<Self::Value, R>
    where R: SerdeError
    {
        let mut encoding = <E::G1Affine as GroupEncoding>::Repr::default();
        encoding.as_mut().copy_from_slice(v);
        let elem = E::G1Affine::from_bytes(&encoding);
        if elem.is_none().into() {
            Err(R::invalid_value(Unexpected::Bytes(v), &self))
        } else {
            Ok(elem.unwrap())
        }
    }
    
}

#[cfg(feature = "serde_support")]
impl<'de, E: Engine> Deserialize<'de> for G1Affine<E> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>
    {
        let visitor = G1AffineVisitor::<E>(PhantomData);
        let inner = deserializer.deserialize_bytes(visitor)?;
        
        Ok(G1Affine(inner))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd)]
pub struct G2Affine<E: Engine>(E::G2Affine);

impl<E: Engine> Copy for G2Affine<E> {}

impl<E: Engine> AsRef<E::G2Affine> for G2Affine<E> {
    fn as_ref(&self) -> &E::G2Affine {
        &self.0
    }
}

impl<E: Engine> AsMut<E::G2Affine> for G2Affine<E> {
    fn as_mut(&mut self) -> &mut E::G2Affine {
        &mut self.0
    }
}

impl<E: Engine> Deref for G2Affine<E> {
    type Target = E::G2Affine;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<E: Engine> G2Affine<E> {
    pub fn into_inner(self) -> E::G2Affine {
        self.0
    }

    pub fn inner(&self) -> &E::G2Affine {
        &self.0
    }

    pub fn from_inner(inner: E::G2Affine) -> Self {
        Self(inner)
    }
}

#[cfg(feature = "serde_support")]
impl<E: Engine> Serialize for G2Affine<E> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer
    {
        serializer.serialize_bytes(self.0.to_bytes().as_ref())
    }
}

#[cfg(feature = "serde_support")]
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd)]
struct G2AffineVisitor<E>(PhantomData<E>);

#[cfg(feature = "serde_support")]
impl<'de, E: Engine> Visitor<'de> for G2AffineVisitor<E> {
    type Value = E::G2Affine;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        write!(formatter, "canonical bytes representation of an affine G1 Element")
    }

    fn visit_bytes<R>(self, v: &[u8]) -> Result<Self::Value, R>
    where R: SerdeError
    {
        let mut encoding = <E::G2Affine as GroupEncoding>::Repr::default();
        encoding.as_mut().copy_from_slice(v);
        let elem = E::G2Affine::from_bytes(&encoding);
        if elem.is_none().into() {
            Err(R::invalid_value(Unexpected::Bytes(v), &self))
        } else {
            Ok(elem.unwrap())
        }
    }
    
}

#[cfg(feature = "serde_support")]
impl<'de, E: Engine> Deserialize<'de> for G2Affine<E> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>
    {
        let visitor = G2AffineVisitor::<E>(PhantomData);
        let inner = deserializer.deserialize_bytes(visitor)?;
        
        Ok(G2Affine(inner))
    }
}

pub struct Scalar<E: Engine>(E::Fr);


#[cfg(all(test, feature = "serde_support"))]
mod tests {
    use bincode::{serialize, deserialize};
    use pairing::group::Curve;
    use rand::{rngs::SmallRng, SeedableRng, Rng};
    use super::{G1Affine, G2Affine};
    use bls12_381::{Bls12, G1Affine as G, G2Affine as H, Scalar};

    #[test]
    fn test_group_elem_serde_support() {
        let seed = [69; 32];
        let mut rng = SmallRng::from_seed(seed);
        let mut g = G::generator();
        let mut h = H::generator();
        let mut gs = vec![G1Affine::<Bls12>::from_inner(g)];
        let mut hs = vec![G2Affine::<Bls12>::from_inner(h)];
        for _ in 0..16 {
            let multiplier: Scalar = rng.gen::<u64>().into();
            g = (g * multiplier).to_affine();
            h = (h * multiplier).to_affine();
            gs.push(G1Affine::from_inner(g));
            hs.push(G2Affine::from_inner(h));
        }

        for g in gs {
            let serialized = serialize(&g).unwrap();
            let deserialized: G1Affine<Bls12> = deserialize(serialized.as_slice()).unwrap();
            assert_eq!(g.into_inner(), deserialized.into_inner())
        }

        for h in hs {
            let serialized = serialize(&h).unwrap();
            let deserialized: G2Affine<Bls12> = deserialize(serialized.as_slice()).unwrap();
            assert_eq!(h.into_inner(), deserialized.into_inner())
        }
    }

}