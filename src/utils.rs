use ark_crypto_primitives::commitment::pedersen::{Commitment, Parameters, Randomness, Window};
use ark_crypto_primitives::commitment::CommitmentScheme;
use ark_ec::CurveGroup;
use ark_ff::{BigInteger, Field, PrimeField};

pub trait ScalarToBytes<C> {
    fn to_bytes(&self) -> Vec<u8>;
}

impl<C: CurveGroup> ScalarToBytes<C> for C::ScalarField {
    fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = vec![];
        for fp in self.to_base_prime_field_elements() {
            bytes.extend_from_slice(&fp.into_bigint().to_bytes_be());
        }
        bytes
    }
}

pub trait CommitVec<C: CurveGroup> {
    fn commit_vec(params: &Parameters<C>, vec: Vec<C::ScalarField>, r: &Randomness<C>)
        -> C::Affine;
}

impl<C: CurveGroup, W: Window> CommitVec<C> for Commitment<C, W> {
    fn commit_vec(
        params: &Parameters<C>,
        vec: Vec<<C>::ScalarField>,
        r: &Randomness<C>,
    ) -> C::Affine {
        let mut bytes = vec![];

        for scalar in vec.iter() {
            bytes.extend_from_slice(&<C::ScalarField as ScalarToBytes<C>>::to_bytes(scalar));
        }

        Self::commit(params, &bytes, r).unwrap()
    }
}
