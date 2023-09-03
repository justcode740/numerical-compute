use nalgebra::{DMatrix, OVector, Dynamic, OMatrix};

pub struct SVDResult {
    pub u: DMatrix<f64>,
    pub singular_values: OVector<f64, Dynamic>,
    pub v_t: DMatrix<f64>,
}

pub fn perform_svd(matrix: OMatrix<f64, Dynamic, Dynamic>) -> SVDResult {
    let svd = matrix.svd(true, true);
    SVDResult {
        u: svd.u.unwrap(),
        singular_values: svd.singular_values,
        v_t: svd.v_t.unwrap(),
    }
}

#[test]
fn test_svd() {
   
    let matrix = OMatrix::<f64, Dynamic, Dynamic>::from_vec(3,3,vec![
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0]
    );

    let svd_result = perform_svd(matrix);

    assert_eq!(svd_result.singular_values[0], 1.0683695145547085_f64);
    assert_eq!(svd_result.singular_values[1], 16.848103352614213_f64);
    assert_eq!(svd_result.singular_values[2], 0_f64);

    assert_eq!(svd_result.u.nrows(), 3);
    assert_eq!(svd_result.singular_values.len(), 3);
    assert_eq!(svd_result.v_t.ncols(), 3);
}
