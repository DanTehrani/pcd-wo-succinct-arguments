use crate::CommitVec;
use ark_crypto_primitives::commitment::pedersen::{Commitment, Parameters, Randomness, Window};
use ark_crypto_primitives::commitment::CommitmentScheme;
use ark_ec::CurveGroup;
use ark_ff::{One, Zero};
use ark_std::{test_rng, UniformRand};
use std::marker::PhantomData;

pub struct HadamardPredicate<C: CurveGroup, W: Window> {
    pub params: Parameters<C>,
    _marker: PhantomData<W>,
}

pub struct HadamardProof<C: CurveGroup> {
    pub c1: C::Affine,
    pub c2: C::Affine,
    pub c3: C::Affine,
    pub w1: C::ScalarField,
    pub w2: C::ScalarField,
    pub w3: C::ScalarField,
}

impl<C: CurveGroup, W: Window> HadamardPredicate<C, W> {
    pub fn setup() -> Self {
        let mut rng = test_rng();
        let params = Commitment::<C, W>::setup(&mut rng).unwrap();
        HadamardPredicate {
            params,
            _marker: PhantomData,
        }
    }

    pub fn hadamard_prod(a: &[C::ScalarField], b: &[C::ScalarField]) -> Vec<C::ScalarField> {
        let mut c = vec![C::ScalarField::zero(); a.len()];
        for i in 0..a.len() {
            c[i] = a[i] * b[i];
        }
        c
    }

    pub fn prove(&self, a: Vec<C::ScalarField>, b: Vec<C::ScalarField>) -> HadamardProof<C> {
        let c = Self::hadamard_prod(&a, &b);
        let mut rng = test_rng();

        let w1 = Randomness::<C>::rand(&mut rng);
        let w2 = Randomness::<C>::rand(&mut rng);
        let w3 = Randomness::<C>::rand(&mut rng);

        let c1 = Commitment::<C, W>::commit_vec(&self.params, a, &w1);
        let c2 = Commitment::<C, W>::commit_vec(&self.params, b, &w2);
        let c3 = Commitment::<C, W>::commit_vec(&self.params, c, &w3);

        HadamardProof {
            c1,
            c2,
            c3,
            w1: w1.0,
            w2: w2.0,
            w3: w3.0,
        }
    }
}
