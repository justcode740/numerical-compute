pub mod svd;
// pub mod eigen_decompose;
// pub mod qr_factorization;

fn main() {
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
}
