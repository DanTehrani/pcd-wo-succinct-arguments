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
use ark_ff::{BigInteger, Field, One, PrimeField, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_std::test_rng;
use std::marker::PhantomData;

pub struct AccumulationProver<C, W>
where
    C: CurveGroup,
    C::BaseField: PrimeField,
{
    pub params: Parameters<C>,
    sponge: PoseidonSponge<C::BaseField>,
    _marker: PhantomData<W>,
}

impl<C, W> AccumulationProver<C, W>
where
    C: CurveGroup,
    C::BaseField: PrimeField,
    W: Window,
{
    pub fn new(sponge: PoseidonSponge<C::BaseField>) -> Self {
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
    ) -> (HadamardAcc<C>, HadamardAccProof<C>) {
        let n = qx.len();
        let l = qw[0].a_vec.len();

        // Absorb the accumulator instances
        for qx_i in &qx {
            self.sponge.absorb(&qx_i.to_absorbable_bytes());
        }

        let mu = self.sponge.squeeze_field_elements::<C::ScalarField>(1)[0];
        let mut mu_powers = vec![];
        for i in 0..n {
            mu_powers.push(mu.pow(&[i as u64]));
        }

        let mut t_vecs = vec![Vec::with_capacity(l); 2 * n - 1];
        for i in 0..l {
            let mut a_coeffs = vec![];
            let mut b_coeffs = vec![];
            for (j, qw_i) in qw.iter().enumerate() {
                a_coeffs.push(qw_i.a_vec[i] * &mu_powers[j]);
                b_coeffs.push(qw_i.b_vec[i]);
            }
            b_coeffs.reverse();

            let a_poly = DensePolynomial::from_coefficients_vec(a_coeffs);
            let b_poly = DensePolynomial::from_coefficients_vec(b_coeffs);

            let product_poly = a_poly.naive_mul(&b_poly);

            let mut product_coeffs = product_poly.coeffs;
            if product_coeffs.len() < 2 * n - 1 {
                product_coeffs.resize_with(2 * n - 1, || C::ScalarField::zero());
            }

            for i in 0..(2 * n - 1) {
                t_vecs[i].push(product_coeffs[i].clone());
                // t_vecs[i]: push the i'th degree coefficient of the product polynomial to the i'th t_vec
                // each t_vec in t_vecs will be full with the coefficients of the product polynomial
                // the first product polynomial will be the degree-0 coefficient of the t_vec polynomials
                // the n'th product polynomial will be the degree-n coefficient of the t_vec polynomials
            }
        }

        // Commit t_vecs
        let mut comm_t_vecs_low = Vec::with_capacity(n - 1);
        let mut comm_t_vecs_high = Vec::with_capacity(n - 1);
        for (i, t_vec) in t_vecs.iter().enumerate() {
            if i == n - 1 {
                continue;
            }

            if i < n {
                let c_t_i = Commitment::<C, W>::commit_vec(
                    &self.params,
                    t_vec.clone(),
                    &Randomness(C::ScalarField::zero()),
                );
                comm_t_vecs_low.push(c_t_i);
            } else {
                let c_t_i = Commitment::<C, W>::commit_vec(
                    &self.params,
                    t_vec.clone(),
                    &Randomness(C::ScalarField::zero()),
                );
                comm_t_vecs_high.push(c_t_i);
            }
        }

        let nu = self.sponge.squeeze_field_elements::<C::ScalarField>(1)[0];

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
        for i in 0..comm_t_vecs_low.len() {
            c3_1 += comm_t_vecs_low[i] * nu_powers[i];
        }
        let mut c3_2 = C::zero();
        for (i, qx_i) in qx.iter().enumerate() {
            c3_2 += qx_i.2 * mu_powers[i];
        }
        c3_2 *= nu_powers[n - 1];

        let mut c3_3 = C::zero();
        for i in 0..comm_t_vecs_high.len() {
            c3_3 += comm_t_vecs_high[i] * nu.pow([(n + i) as u64]);
        }

        let c3 = c3_1 + c3_2 + c3_3;

        // a_1 * mu^0 * nu^0 + a_2 * mu^1 * nu^1
        let mut a = vec![C::ScalarField::zero(); l];
        let mut b = vec![C::ScalarField::zero(); l];

        // qw_1.w_1 * mu^0 * nu^0 + qw_1.w_2 * mu^1 * nu^1
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

        let mut comm_t_vecs = vec![];
        comm_t_vecs.extend_from_slice(&comm_t_vecs_low);
        comm_t_vecs.extend_from_slice(&comm_t_vecs_high);

        (
            HadamardAcc {
                qx: HadamardInstance(c1.into_affine(), c2.into_affine(), c3.into_affine()),
                qw: HadamardWitness {
                    a_vec: a,
                    b_vec: b,
                    w1,
                    w2,
                    w3,
                    _marker: PhantomData,
                },
            },
            HadamardAccProof(comm_t_vecs),
        )
    }
}
