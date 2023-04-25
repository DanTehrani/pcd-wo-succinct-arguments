use crate::hadamard::accumulation::{HadamardAccProof, HadamardInstance};
use ark_crypto_primitives::commitment::pedersen::Window;
use ark_crypto_primitives::sponge::{
    poseidon::PoseidonSponge, CryptographicSponge, FieldBasedCryptographicSponge,
};
use ark_ec::{CurveGroup, Group};
use ark_ff::{Field, PrimeField};
use std::marker::PhantomData;

pub struct AccumulationVerifier<C, W>
where
    C: CurveGroup,
    C::BaseField: PrimeField,
    W: Window,
{
    sponge: PoseidonSponge<C::BaseField>,
    _marker: PhantomData<W>,
}

impl<C, W> AccumulationVerifier<C, W>
where
    C: CurveGroup,
    C::BaseField: PrimeField,
    W: Window,
{
    pub fn new(sponge: PoseidonSponge<C::BaseField>) -> Self {
        AccumulationVerifier {
            sponge,
            _marker: PhantomData,
        }
    }

    pub fn verify(
        &mut self,
        acc_x: &HadamardInstance<C>,
        acc_xs: &Vec<HadamardInstance<C>>,
        proof: &HadamardAccProof<C>,
    ) {
        // Absorb the accumulator instances
        for acc_xi in acc_xs {
            self.sponge.absorb(&acc_xi.to_absorbable_bytes());
        }

        let n = acc_xs.len();
        let mu: C::ScalarField = self.sponge.squeeze_field_elements(1)[0];
        let nu: C::ScalarField = self.sponge.squeeze_field_elements(1)[0];

        let mut mu_powers = vec![];
        for i in 0..n {
            mu_powers.push(mu.pow(&[i as u64]));
        }

        let mut nu_powers = vec![];
        for i in 0..n {
            nu_powers.push(nu.pow(&[i as u64]));
        }

        let mut expected_c1 = C::zero();
        let mut expected_c2 = C::zero();
        for i in 0..n {
            expected_c1 += acc_xs[i].0 * mu_powers[i] * nu_powers[i];
            expected_c2 += acc_xs[i].1 * nu_powers[nu_powers.len() - 1 - i];
        }

        let mut expected_c3_1 = C::zero();
        let mut expected_c3_2 = C::zero();
        let mut expected_c3_3 = C::zero();

        for i in 0..(n - 1) {
            expected_c3_1 += proof.0[i] * nu_powers[i];
        }

        for i in 0..acc_xs.len() {
            expected_c3_2 += acc_xs[i].2 * mu_powers[i];
        }
        expected_c3_2 *= nu_powers[n - 1];

        for i in 0..(n - 1) {
            expected_c3_3 += proof.0[i + n - 1] * nu.pow([(n + i) as u64]);
        }

        let expected_c3 = expected_c3_1 + expected_c3_2 + expected_c3_3;

        assert_eq!(expected_c1.into_affine(), acc_x.0);
        assert_eq!(expected_c2.into_affine(), acc_x.1);
        assert_eq!(expected_c3.into_affine(), acc_x.2);
    }
}
