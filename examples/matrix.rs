use std::env::{var, set_var};
use clap::Parser;
use halo2_scaffold::scaffold::cmd::Cli;
use halo2_scaffold::scaffold::run_builder_on_inputs;
use nalgebra::{ Matrix2};

use halo2_base::utils::{ScalarField};
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
    pub u: Vec<Vec<f64>>, 
    pub s: Vec<Vec<f64>>,
    pub vt: Vec<Vec<f64>>,
    pub matrix: Vec<Vec<f64>>,
    pub err: f64 // err tolerance for verification
}

// this algorithm takes a public input x, computes x^2 + 72, and outputs the result as public output
fn some_algorithm_in_zk<F: ScalarField>(
    ctx: &mut Context<F>,
    input: CircuitInput,
    make_public: &mut Vec<AssignedValue<F>>,
) {
    

    let lookup_bits =
        var("LOOKUP_BITS").unwrap_or_else(|_| panic!("LOOKUP_BITS not set")).parse().unwrap();
    const PRECISION_BITS: u32 = 32;
    // fixed-point exp arithmetic
    let fixed_point_chip = FixedPointChip::<F, PRECISION_BITS>::default(lookup_bits);

    let res = naive_matrix_mul(&naive_matrix_mul(&input.u, &input.s), &input.vt);
    res.iter().for_each(|row| {
        println!("res: {}", row.iter().map(|&e| e.to_string()).collect::<Vec<_>>().join(" "));
    });
    
    let mut acc = ctx.load_zero();
    for i in 0..res.len() {
        for j in 0..res[0].len() {
            let val = fixed_point_chip.quantization(res[i][j]);
            let target = fixed_point_chip.quantization(input.matrix[i][j]);
            let val = ctx.load_witness(val);
            let target = ctx.load_witness(target);
            let diff = fixed_point_chip.qsub(ctx, val, target);
            print!("diff val: {:?}", fixed_point_chip.dequantization(*diff.value()));
            let diff_square = fixed_point_chip.qmul(ctx, diff, diff);
            print!("diff sqaure val: {:?}", fixed_point_chip.dequantization(*diff_square.value()));
            acc = fixed_point_chip.qadd(ctx, acc, diff_square);
            print!("acc val: {:?}", fixed_point_chip.dequantization(*acc.value()));
        }
    }
    
    let err = fixed_point_chip.qsqrt(ctx, acc);
    let err_decimal = fixed_point_chip.dequantization(*err.value());
    println!("err: {:?}", err_decimal);
    assert!(err_decimal < input.err); 
   
}

fn naive_matrix_mul(u: &Vec<Vec<f64>>, s: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let num_rows_u = u.len();
    let num_cols_u = u[0].len();
    let num_rows_s = s.len();
    let num_cols_s = s[0].len();

    assert_eq!(num_cols_u, num_rows_s, "Dimensions mismatch for multiplication.");

    let mut result = vec![vec![0.0; num_cols_s]; num_rows_u];

    for i in 0..num_rows_u {
        for j in 0..num_cols_s {
            for k in 0..num_cols_u {
                result[i][j] += u[i][k] * s[k][j];
            }
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
    
    let matrix2 = Matrix2::from_vec(
        matrix.concat()
    );

    // Perform SVD
    let svd = matrix2.svd(true, true);
  
    // Convert U to Vec<Vec<f64>>
    let u_vec: Vec<Vec<f64>> = vec![
        vec![svd.u.unwrap().m11, svd.u.unwrap().m12],
        vec![svd.u.unwrap().m21, svd.u.unwrap().m22]
    ];

    let s = vec![svd.singular_values[0], svd.singular_values[1]];
    let mut eigen_matrix: Vec<Vec<f64>> = vec![vec![0.0; 2]; 2];

    // Fill the diagonal elements with values from the vector
    for i in 0..2 {
        eigen_matrix[i][i] = s[i];
    }

    // Convert V^T to Vec<Vec<f64>>
    let v_t_vec: Vec<Vec<f64>> = vec![
        vec![svd.v_t.unwrap().m11, svd.v_t.unwrap().m12],
        vec![svd.v_t.unwrap().m21, svd.v_t.unwrap().m22]
    ];

    u_vec.iter().for_each(|row| {
        println!("uvec: {}", row.iter().map(|&e| e.to_string()).collect::<Vec<_>>().join(" "));
    });

    eigen_matrix.iter().for_each(|row| {
        println!("eigen: {}", row.iter().map(|&e| e.to_string()).collect::<Vec<_>>().join(" "));
    });

    v_t_vec.iter().for_each(|row| {
        println!("vt: {}", row.iter().map(|&e| e.to_string()).collect::<Vec<_>>().join(" "));
    });

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

    let input = CircuitInput { u: m1, s: m2, vt: m3, matrix, err: 0.5  };
    // run different zk commands based on the command line arguments
    run_builder_on_inputs(
        |builder, input, public| some_algorithm_in_zk(builder.main(0), input, public),
        args,
        input,
    );
}
