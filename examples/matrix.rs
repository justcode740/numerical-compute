use std::env::{var, set_var};
use axiom_eth::{EthChip, keccak::KeccakChip, Field};
use clap::Parser;
use halo2_scaffold::scaffold::{cmd::Cli, run_eth_builder_on_inputs};
use halo2_scaffold::scaffold::{run_builder_on_inputs, run_eth};
use halo2_base::{utils::ScalarField, gates::{GateChip, GateInstructions}};
use poseidon::PoseidonChip;

use halo2_base::AssignedValue;
#[allow(unused_imports)]
use halo2_base::{
    Context,
    QuantumCell::{Constant, Existing, Witness},
};
use serde::{Deserialize, Serialize};
use zkfixedpointchip::gadget::fixed_point::{FixedPointChip, FixedPointInstructions};


#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CircuitInput {
    pub u: Vec<Vec<f64>>, //m * m
    pub s: Vec<Vec<f64>>, //m*n
    pub vt: Vec<Vec<f64>>,//n*n
    pub matrix: Vec<Vec<f64>>, // m*n
    pub err: f64 // err tolerance for verification
}

fn verify_svd_random<F: Field>(
    ctx: &mut Context<F>,
    eth_chip: &EthChip<F>,
    keccak: &mut KeccakChip<F>,
    input: CircuitInput,
    make_public: &mut Vec<AssignedValue<F>>,
) -> impl FnOnce(&mut Context<F>, &mut Context<F>, &EthChip<F>) + Clone {
    let m = input.matrix.len();
    let n = input.matrix[0].len();

    let lookup_bits =
        var("LOOKUP_BITS").unwrap_or_else(|_| panic!("LOOKUP_BITS not set")).parse().unwrap();
    const PRECISION_BITS: u32 = 32;
    // fixed-point exp arithmetic
    let fixed_point_chip = FixedPointChip::<F, PRECISION_BITS>::default(lookup_bits);

    let mut assgined_matrix: Vec<Vec<AssignedValue<F>>> = Vec::new();
    let mut assgined_u: Vec<Vec<AssignedValue<F>>> = Vec::new();
    let mut assgined_s: Vec<Vec<AssignedValue<F>>> = Vec::new();
    let mut assgined_vt: Vec<Vec<AssignedValue<F>>> = Vec::new();

    for i in 0..m {
        let mut row1 = Vec::new(); 
        let mut row2 = Vec::new();
        for j in 0..n {
            let val = ctx.load_witness(fixed_point_chip.quantization(input.matrix[i][j]));
            let eigen = ctx.load_witness(fixed_point_chip.quantization(input.s[i][j]));
            make_public.push(val);
            make_public.push(eigen);
            row1.push(val);
            row2.push(eigen);
        }
        assgined_matrix.push(row1); 
        assgined_s.push(row2);
    }

    for i in 0..m {
        let mut row: Vec<AssignedValue<F>> = Vec::new(); 
        for j in 0..m {
            let val = ctx.load_witness(fixed_point_chip.quantization(input.u[i][j]));
            make_public.push(val);
            row.push(val);
        }
        assgined_u.push(row); 
    }

    for i in 0..n {
        let mut row: Vec<AssignedValue<F>> = Vec::new(); 
        for j in 0..n {
            let val = ctx.load_witness(fixed_point_chip.quantization(input.vt[i][j]));
            make_public.push(val);
            row.push(val);
        }
        assgined_vt.push(row); 
    }

    const T: usize = 3;
    const RATE: usize = 2;
    const R_F: usize = 8;
    const R_P: usize = 57;
    #[allow(clippy::let_and_return)]
    let callback =
        move |ctx_gate: &mut Context<F>, ctx_rlc: &mut Context<F>, eth_chip: &EthChip<F>| {
            let random = eth_chip.rlc().gamma();
            println!("r value is {:?}", random);
            let mut random_ass = ctx_rlc.load_witness(*random);
            let gate = GateChip::<F>::default();
            let mut random_vec = vec![];
            let mut poseidon : PoseidonChip<F, T, RATE> = PoseidonChip::<F, T, RATE>::new(ctx_rlc, R_F, R_P).unwrap();
            // poseidon = PoseidonChip::new(ctx_gate, spec, range_chip);
            // construct random vec n * 1
            // for _ in 0..n/32 {
            //     let bits : Vec<AssignedValue<F>> = gate.num_to_bits(ctx_gate, random_ass, 32);
            //     random_vec.extend(bits);
                
            //     poseidon.update(&[random_ass]);
            //     random_ass = poseidon.squeeze(ctx_gate, &gate).unwrap();
            // }
            // println!("{:?} {:?}", random_vec.len(), 32 * (n/32));
            // assert!(random_vec.len() == 32 * (n/32), "len of random vector incorrect");
    
            // let bits = eth_chip.gate().num_to_bits(ctx_rlc, random_ass, 32);
            for i in 0..n%32 {
                random_vec.push(random_ass);
            }
            assert!(random_vec.len() == n, "len of random vector incorrect");

            // compute u s vt apply to v and matrix apply to v
            let vtv = naive_matrix_vec_mul(ctx_rlc, &fixed_point_chip, &assgined_vt, &random_vec);
            let svtv = naive_matrix_vec_mul(ctx_rlc, &fixed_point_chip, &assgined_s, &vtv);
            let usvtv = naive_matrix_vec_mul(ctx_rlc, &fixed_point_chip, &assgined_u, &svtv);
            
            // let matrixv = naive_matrix_vec_mul(ctx_rlc, &fixed_point_chip, &assgined_matrix, &random_vec);


            // // compute err
            // let mut acc = ctx_rlc.load_zero();
            // for i in 0..m {
            //     let diff = fixed_point_chip.qsub(ctx_rlc, usvtv[i], matrixv[i]);
            //     let abs_diff = fixed_point_chip.qabs(ctx_rlc, diff);
            //     acc = fixed_point_chip.qadd(ctx_rlc, acc, abs_diff);
            // }
            // let matrix_size = ctx_rlc.load_constant(fixed_point_chip.quantization((m * n) as f64));
            // let err = fixed_point_chip.qdiv(ctx_rlc, acc, matrix_size);
            // let err_decimal = fixed_point_chip.dequantization(*err.value());
            // println!("err: {:?}", err_decimal);
            // assert!(err_decimal < input.err);
        };
        
    callback 
   
}

// this algorithm takes a public input 
fn verify_svd_naive<F: ScalarField>(
    ctx: &mut Context<F>,
    input: CircuitInput,
    make_public: &mut Vec<AssignedValue<F>>,
) {
    let m = input.matrix.len();
    let n = input.matrix[0].len();

    let lookup_bits =
        var("LOOKUP_BITS").unwrap_or_else(|_| panic!("LOOKUP_BITS not set")).parse().unwrap();
    const PRECISION_BITS: u32 = 32;
    // fixed-point exp arithmetic
    let fixed_point_chip = FixedPointChip::<F, PRECISION_BITS>::default(lookup_bits);

    let mut assgined_matrix: Vec<Vec<AssignedValue<F>>> = Vec::new();
    let mut assgined_u: Vec<Vec<AssignedValue<F>>> = Vec::new();
    let mut assgined_s: Vec<Vec<AssignedValue<F>>> = Vec::new();
    let mut assgined_vt: Vec<Vec<AssignedValue<F>>> = Vec::new();

    for i in 0..m {
        let mut row1 = Vec::new(); 
        let mut row2 = Vec::new();
        for j in 0..n {
            let val = ctx.load_witness(fixed_point_chip.quantization(input.matrix[i][j]));
            let eigen = ctx.load_witness(fixed_point_chip.quantization(input.s[i][j]));
            make_public.push(val);
            make_public.push(eigen);
            row1.push(val);
            row2.push(eigen);
        }
        assgined_matrix.push(row1); 
        assgined_s.push(row2);
    }

    for i in 0..m {
        let mut row: Vec<AssignedValue<F>> = Vec::new(); 
        for j in 0..m {
            let val = ctx.load_witness(fixed_point_chip.quantization(input.u[i][j]));
            make_public.push(val);
            row.push(val);
        }
        assgined_u.push(row); 
    }

    for i in 0..n {
        let mut row: Vec<AssignedValue<F>> = Vec::new(); 
        for j in 0..n {
            let val = ctx.load_witness(fixed_point_chip.quantization(input.vt[i][j]));
            make_public.push(val);
            row.push(val);
        }
        assgined_vt.push(row); 
    }

    let us = naive_matrix_mul(ctx, & fixed_point_chip, &assgined_u, &assgined_s);
    let usvt = naive_matrix_mul(ctx, & fixed_point_chip, &us, &assgined_vt);

    let mut acc = ctx.load_zero();
    for i in 0..m {
        for j in 0..n {
            let diff = fixed_point_chip.qsub(ctx, usvt[i][j], assgined_matrix[i][j]);
            let abs_diff = fixed_point_chip.qabs(ctx, diff);
            acc = fixed_point_chip.qadd(ctx, acc, abs_diff);
        }
    }
    let matrix_size = ctx.load_constant(fixed_point_chip.quantization((m * n) as f64));
    let err = fixed_point_chip.qdiv(ctx, acc, matrix_size);
    let err_decimal = fixed_point_chip.dequantization(*err.value());
    println!("err: {:?}", err_decimal);
    assert!(err_decimal < input.err); 
   
}

fn naive_matrix_mul<F: ScalarField>(ctx: &mut Context<F>, fixed_point_chip: &FixedPointChip<F, 32>, u: &Vec<Vec<AssignedValue<F>>>, s: &Vec<Vec<AssignedValue<F>>>) ->  Vec<Vec<AssignedValue<F>>> {
    let num_rows_u = u.len();
    let num_cols_u = u[0].len();
    let num_rows_s = s.len();
    let num_cols_s = s[0].len();

    assert_eq!(num_cols_u, num_rows_s, "Dimensions mismatch for multiplication.");

    let mut result = vec![vec![ctx.load_zero(); num_cols_s]; num_rows_u];

    for i in 0..num_rows_u {
        for j in 0..num_cols_s {
            let mut acc = ctx.load_constant(F::from(0));
            for k in 0..num_cols_u {
                let mul_res = fixed_point_chip.qmul(ctx, u[i][k], s[k][j]);
                acc = fixed_point_chip.qadd(ctx, mul_res, acc);
            }
            result[i][j] = acc;
        }
    }

    result
}

fn naive_matrix_vec_mul<F: ScalarField>(ctx: &mut Context<F>, fixed_point_chip: &FixedPointChip<F, 32>, matrix: &Vec<Vec<AssignedValue<F>>>,vec: &Vec<AssignedValue<F>>) ->  Vec<AssignedValue<F>> {
    let n = vec.len();
    let m = matrix.len();
    assert_eq!(matrix[0].len(), n, "Dimensions mismatch for multiplication.");
    let mut result = vec![ctx.load_zero(); m];
    for i in 0..m {
        let mut acc = ctx.load_constant(F::from(0));
        for j in 0..n {
            let interm = fixed_point_chip.qmul(ctx, matrix[i][j], vec[j]);
            acc = fixed_point_chip.qadd(ctx, acc, interm);
        }
        result[i] = acc;
    }

    result
}

// assume square matrix, later padding with zero
fn strassen_matrix_mul<F: ScalarField>(ctx: &mut Context<F>, fixed_point_chip: &FixedPointChip<F, 32> ,u: &Vec<Vec<AssignedValue<F>>>, s: &Vec<Vec<AssignedValue<F>>>) ->  Vec<Vec<AssignedValue<F>>> {
    // ctx.constrain_equal(a, b)
    let num_rows_u: usize = u.len();
    let num_cols_u = u[0].len();
    let num_rows_s = s.len();
    let num_cols_s = s[0].len();

    if num_rows_s > num_cols_u {
        
    }

    assert_eq!(num_cols_u, num_rows_s, "Dimensions mismatch for multiplication.");

    let mut result = vec![vec![ctx.load_zero(); num_cols_s]; num_rows_u];

    for i in 0..num_rows_u {
        for j in 0..num_cols_s {
            let mut acc = ctx.load_constant(F::from(0));
            for k in 0..num_cols_u {
                let mul_res = fixed_point_chip.qmul(ctx, u[i][k], s[k][j]);
                acc = fixed_point_chip.qadd(ctx, mul_res, acc);
            }
            result[i][j] = acc;
        }
    }

    result
}



fn main() {
    env_logger::init();
    set_var("LOOKUP_BITS", 12.to_string());
    set_var("DEGREE", 13.to_string());

    let matrix = vec![
        vec![1.0,2.0],
        vec![3.0,4.0]
    ];
    
    let args = Cli::parse();

    let m1: Vec<Vec<f64>> = vec![
        vec![-0.40455358, -0.9145143],
        vec![-0.9145143, 0.40455358],
    ];

    let m2: Vec<Vec<f64>> = vec![
        vec![5.4649857, 0.0],
        vec![0.0, 0.36596619]
    ];

    let m3: Vec<Vec<f64>> = vec![
        vec![-0.57604844, -0.81741556],
        vec![0.81741556, -0.57604844],
    ];

    let input = CircuitInput { u: m1, s: m2, vt: m3, matrix, err: 1e-5 };
    // run different zk commands based on the command line arguments
    // run_builder_on_inputs(
    //     |builder, input, public| verify_svd_naive(builder.main(0), input, public),
    //     args,
    //     input
    // );
    // run_eth(|ctx, ethchip, keccakchip, input, public| verify_svd_random(ctx, ethchip, keccakchip, input, public), args);
    run_eth_builder_on_inputs(|builder, ethchip, keccakchip, input, public| verify_svd_random(builder.main(0), ethchip, keccakchip, input, public),
    args,
    input
    )
}
