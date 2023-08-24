# Numerical Linear Algebra in Halo2

the goal is to implement proof system based on halo2 to verify linear algebra operations
like svd, qr, eigen decomposition.

### setup

`src/reference`: implementation of algebra computation in plain rust

`src/zk`: halo2 zk circuits

run `cargo test` in root folder to run tests
run `cargo test` in `src/zk` for zk specific tests

### todo 
- impl svd in halo2 without too much optimization
    - how to do fp
    - unconstrained function?
- testing
- e2e proof system, run some metrics
- optimization
- add qr and eigen decomposition
- how to do fs with halo2 for randomness check protocol: https://hackmd.io/@axiom/SJw3p-qX3





