mod decider;
mod prover;
mod verifier;

use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{BigInteger, Field, PrimeField};
use std::marker::PhantomData;

pub struct HadamardAccProof<C>(Vec<C::Affine>)
where
    C: CurveGroup;

#[derive(Clone)]
pub struct HadamardInstance<C: CurveGroup>(C::Affine, C::Affine, C::Affine);

impl<C: CurveGroup> HadamardInstance<C> {
    fn to_absorbable_bytes(&self) -> Vec<u8> {
        let c_1 = self.0;
        let c_2 = self.1;
        let c_3 = self.2;
        let (c1_x, c1_y) = c_1.xy().unwrap();
        let (c2_x, c2_y) = c_2.xy().unwrap();
        let (c3_x, c3_y) = c_3.xy().unwrap();

        let mut bytes = vec![];
        for field_element in [c1_x, c1_y, c2_x, c2_y, c3_x, c3_y].iter() {
            for fp in field_element.to_base_prime_field_elements() {
                bytes.extend_from_slice(&fp.into_bigint().to_bytes_be());
            }
        }

        bytes
    }
}

#[derive(Clone)]
pub struct HadamardWitness<C>
where
    C: CurveGroup,
{
    pub a_vec: Vec<C::ScalarField>,
    pub b_vec: Vec<C::ScalarField>,
    pub w1: C::ScalarField,
    pub w2: C::ScalarField,
    pub w3: C::ScalarField,
    _marker: PhantomData<C>,
}

pub struct HadamardAcc<C>
where
    C: CurveGroup,
{
    pub qx: HadamardInstance<C>,
    pub qw: HadamardWitness<C>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hadamard::accumulation::decider::Decider;
    use crate::hadamard::accumulation::prover::AccumulationProver;
    use crate::hadamard::accumulation::verifier::AccumulationVerifier;
    use crate::hadamard::predicate::{HadamardPredicate, HadamardProof};
    use crate::parameters::*;
    use ark_bls12_381::{g1::G1Projective, Fr};
    use ark_crypto_primitives::commitment::{
        pedersen::{Commitment, Window},
        CommitmentScheme,
    };
    use ark_crypto_primitives::sponge::{
        poseidon::PoseidonSponge, CryptographicSponge, FieldBasedCryptographicSponge,
    };

    #[derive(Clone)]
    struct HadamardCommWindow {}

    impl Window for HadamardCommWindow {
        const WINDOW_SIZE: usize = 256;
        const NUM_WINDOWS: usize = 2;
    }

    #[derive(Clone)]
    struct AccProverCommWindow {}

    impl Window for AccProverCommWindow {
        const WINDOW_SIZE: usize = 256;
        const NUM_WINDOWS: usize = 4;
    }

    #[test]
    // Split accumulate two predicate instances
    fn test_hs_acc() {
        // Setup
        let hadamard_predicate = HadamardPredicate::<G1Projective, HadamardCommWindow>::setup();
        let poseidon_config = get_poseidon_params::<G1Projective>(2);

        let a_1 = vec![Fr::from(1), Fr::from(2)];
        let b_1 = vec![Fr::from(2), Fr::from(3)];

        let a_2 = vec![Fr::from(1), Fr::from(2)];
        let b_2 = vec![Fr::from(4), Fr::from(5)];

        // Accumulate two predicate instances

        let proof_1 = hadamard_predicate.prove(a_1.clone(), b_1.clone());
        let proof_2 = hadamard_predicate.prove(a_2.clone(), b_2.clone());

        let qx_1 = HadamardInstance::<G1Projective>(proof_1.c1, proof_1.c2, proof_1.c3);
        let qw_1 = HadamardWitness::<G1Projective> {
            a_vec: a_1,
            b_vec: b_1,
            w1: proof_1.w1,
            w2: proof_1.w2,
            w3: proof_1.w3,
            _marker: PhantomData,
        };

        let qx_2 = HadamardInstance::<G1Projective>(proof_2.c1, proof_2.c2, proof_2.c3);
        let qw_2 = HadamardWitness::<G1Projective> {
            a_vec: a_2,
            b_vec: b_2,
            w1: proof_2.w1,
            w2: proof_2.w2,
            w3: proof_2.w3,
            _marker: PhantomData,
        };

        let qx = vec![qx_1, qx_2.clone(), qx_2];
        let qw = vec![qw_1, qw_2.clone(), qw_2];

        // Compute the accumulator and the accumulation proof
        let prover_sponge = PoseidonSponge::new(&poseidon_config);
        let mut prover =
            AccumulationProver::<G1Projective, AccProverCommWindow>::new(prover_sponge);

        let (acc, pf) = prover.prove_acc(qx.clone(), qw);

        let verifier_sponge = PoseidonSponge::new(&poseidon_config);
        let mut verifier =
            AccumulationVerifier::<G1Projective, AccProverCommWindow>::new(verifier_sponge);

        // Verify the accumulation
        verifier.verify(&acc.qx, &qx, &pf);

        let decider = Decider::<G1Projective, HadamardCommWindow>::new(hadamard_predicate.params);
        decider.verify(&acc);
    }
}
