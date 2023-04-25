use crate::hadamard::predicate::HadamardPredicate;
use crate::{hadamard::accumulation::HadamardAcc, CommitVec};
use ark_crypto_primitives::{
    commitment::pedersen::{Commitment, Parameters, Randomness},
    crh::pedersen::Window,
};
use ark_ec::CurveGroup;
use std::marker::PhantomData;

pub struct Decider<C: CurveGroup, W: Window> {
    pub params: Parameters<C>,
    _marker: PhantomData<W>,
}

impl<C: CurveGroup, W: Window> Decider<C, W> {
    pub fn new(params: Parameters<C>) -> Self {
        Self {
            params,
            _marker: PhantomData,
        }
    }

    pub fn verify(&self, acc: &HadamardAcc<C>) {
        let a = acc.qw.clone().a_vec;
        let b = acc.qw.clone().b_vec;

        let w1 = acc.qw.w1;
        let w2 = acc.qw.w2;
        let w3 = acc.qw.w3;

        let c: Vec<C::ScalarField> = HadamardPredicate::<C, W>::hadamard_prod(&a, &b);

        let c_1 = Commitment::<C, W>::commit_vec(&self.params, a, &Randomness(w1));
        let c_2 = Commitment::<C, W>::commit_vec(&self.params, b, &Randomness(w2));
        let c_3 = Commitment::<C, W>::commit_vec(&self.params, c, &Randomness(w3));

        //        assert_eq!(c_1, acc.qx.0);
        ///        assert_eq!(c_2, acc.qx.1);
        assert_eq!(c_3, acc.qx.2);
    }
}
