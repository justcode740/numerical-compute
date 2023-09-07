use std::env::{var, set_var};
use clap::Parser;
use halo2_scaffold::scaffold::cmd::Cli;
use halo2_scaffold::scaffold::run_builder_on_inputs;
use halo2_base::utils::ScalarField;
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

// this algorithm takes a public input x, computes x^2 + 72, and outputs the result as public output
fn some_algorithm_in_zk<F: ScalarField>(
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

fn naive_matrix_mul<F: ScalarField>(ctx: &mut Context<F>, fixed_point_chip: &FixedPointChip<F, 32> ,u: &Vec<Vec<AssignedValue<F>>>, s: &Vec<Vec<AssignedValue<F>>>) ->  Vec<Vec<AssignedValue<F>>> {
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
    run_builder_on_inputs(
        |builder, input, public| some_algorithm_in_zk(builder.main(0), input, public),
        args,
        input,
    );
}
