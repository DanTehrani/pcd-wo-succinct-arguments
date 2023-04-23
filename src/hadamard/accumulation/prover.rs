use crate::hadamard::accumulation::{
    HadamardAcc, HadamardAccProof, HadamardInstance, HadamardWitness,
};
use crate::CommitVec;
use ark_crypto_primitives::commitment::pedersen::{Commitment, Parameters, Randomness, Window};
use ark_crypto_primitives::commitment::CommitmentScheme;
use ark_crypto_primitives::sponge::{
    poseidon::PoseidonSponge, CryptographicSponge, FieldBasedCryptographicSponge,
};
use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_ff::{BigInteger, Field, PrimeField, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_std::test_rng;
use std::marker::PhantomData;

pub struct AccumulationProver<C, W>
where
    C: CurveGroup,
{
    pub params: Parameters<C>,
    sponge: PoseidonSponge<<C as Group>::ScalarField>,
    _marker: PhantomData<W>,
}

impl<C, W> AccumulationProver<C, W>
where
    C: CurveGroup,
    W: Window,
{
    pub fn new(sponge: PoseidonSponge<<C as Group>::ScalarField>) -> Self {
        let mut rng = test_rng();
        let params = Commitment::<C, W>::setup(&mut rng).unwrap();
        AccumulationProver {
            params,
            sponge,
            _marker: PhantomData,
        }
    }

    // Split accumulate Hadamard predicate instances and witnesses
    pub fn prove_acc(
        &mut self,
        qx: Vec<HadamardInstance<C>>,
        qw: Vec<HadamardWitness<C>>,
    ) -> (
        (HadamardInstance<C>, HadamardWitness<C>),
        HadamardAccProof<C>,
    ) {
        let n = qx.len();
        let l = qw[0].a_vec.len();

        // Absorb the accumulator instances
        for qx_i in &qx {
            self.sponge.absorb(&qx_i.to_absorbable_bytes());
        }

        let mu = self.sponge.squeeze_native_field_elements(1)[0];
        let mut mu_powers = vec![];
        for i in 0..n {
            mu_powers.push(mu.pow(&[i as u64]));
        }

        let num_inputs = qx.len();
        let mut t_vecs = vec![Vec::with_capacity(l); 2 * num_inputs - 1];
        for (i, qw_i) in qw.iter().enumerate() {
            let mu_a_vec = qw_i.a_vec.iter().map(|a| *a * &mu_powers[i]).collect();
            let mut b_vec = qw_i.b_vec.clone();
            b_vec.reverse();

            let a_poly = DensePolynomial::from_coefficients_vec(mu_a_vec);
            let b_poly = DensePolynomial::from_coefficients_vec(b_vec);

            let product_poly = a_poly.naive_mul(&b_poly);
            let mut product_coeffs = product_poly.coeffs;
            if product_coeffs.len() < 2 * num_inputs - 1 {
                product_coeffs.resize_with(2 * num_inputs - 1, || C::ScalarField::zero());
            }

            for i in 0..(2 * num_inputs - 1) {
                t_vecs[i].push(product_coeffs[i].clone());
                // t_vecs[i]: push the i'th degree coefficient of the product polynomial to the i'th t_vec
                // each t_vec in t_vecs will be full with the coefficients of the product polynomial
                // the first product polynomial will be the degree-0 coefficient of the t_vec polynomials
                // the n'th product polynomial will be the degree-n coefficient of the t_vec polynomials
            }
        }

        // Commit t_vecs
        let mut comm_t_vecs = Vec::with_capacity(n * 2 - 1);
        for t_vec in t_vecs {
            let c_t_i = Commitment::<C, W>::commit_vec(
                &self.params,
                t_vec,
                &Randomness(C::ScalarField::zero()),
            );
            comm_t_vecs.push(c_t_i);
        }

        let nu = self.sponge.squeeze_native_field_elements(1)[0];

        let mut nu_powers = vec![];
        for i in 0..qx.len() {
            nu_powers.push(nu.pow(&[i as u64]));
        }

        // Compute commitment to a(v, u)
        let mut c1 = C::zero();
        let mut c2 = C::zero();
        for (i, qx_i) in qx.iter().enumerate() {
            c1 += qx_i.0 * nu_powers[i] * mu_powers[i];
            c2 += qx_i.1 * nu_powers[nu_powers.len() - i - 1];
        }

        let mut c3_1 = C::zero();
        for i in 0..n {
            c3_1 += comm_t_vecs[i] * nu_powers[i];
        }
        let mut c3_2 = C::zero();
        for (i, qx_i) in qx.iter().enumerate() {
            c3_2 += qx_i.2 * mu_powers[i];
        }
        c3_2 *= nu_powers[n - 1];

        let mut c3_3 = C::zero();
        for i in 0..(n - 1) {
            c3_3 += comm_t_vecs[n + i] * nu.pow([(n + i - 1) as u64]);
        }

        let c3 = c3_1 + c3_2 + c3_3;

        let mut a = vec![C::ScalarField::zero(); qx.len()];
        let mut b = vec![C::ScalarField::zero(); qx.len()];
        let mut w1 = C::ScalarField::zero();
        let mut w2 = C::ScalarField::zero();
        let mut w3 = C::ScalarField::zero();

        for (i, qw_i) in qw.iter().enumerate() {
            for j in 0..qw_i.a_vec.len() {
                a[j] += qw_i.a_vec[j] * &mu_powers[i] * &nu_powers[i];
            }
            w1 += qw_i.w1 * &mu_powers[i] * &nu_powers[i];

            for j in 0..qw_i.b_vec.len() {
                b[j] += qw_i.b_vec[j] * &nu_powers[n - i - 1];
            }
            w2 += qw_i.w2 * &nu_powers[n - i - 1];

            w3 += qw_i.w3 * &mu_powers[i];
        }
        w3 *= nu_powers[n - 1];

        (
            (
                HadamardInstance(c1.into_affine(), c2.into_affine(), c3.into_affine()),
                HadamardWitness {
                    a_vec: a,
                    b_vec: b,
                    w1,
                    w2,
                    w3,
                    _marker: PhantomData,
                },
            ),
            HadamardAccProof(comm_t_vecs),
        )
    }
}
