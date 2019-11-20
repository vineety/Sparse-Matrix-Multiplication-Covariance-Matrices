# Sparse-Matrix-Multiplication-Covariance-Matrices

Primary Algorithm : Algorithmically, Sparse-Sparse multiplication problems manifests itself in three possible forms:(a) Multiplication of a sparse matrix with a sparse diagonal, sparse block-diagonal, or full-dense covariance matrix ,(b)Multiplication of a sparse matrix with a sparse matrix that result in a sparse symmetric matrix, and (c) Multiplication of a sparse matrix multiplication with a sparse that result in a dense symmetric matrix. The code as part of this repository performs these operations in computationally efficiencint manner.

Problem Formulation:

Matrix multiplication of two sparse matrices is a fundamental operation in linear Bayesian inverse problems for computing covariance matrices of observations and a posteriori uncertainties. Applications of sparse-sparse matrix multiplication algorithms for specific use-cases in such inverse problems remain unexplored. We present a hybrid-parallel sparse-sparse matrix multiplication approach that is more efficient by a third in terms of execution time and operation count relative to standard sparse matrix multiplication algorithms available in most libraries. Two modifications of this hybrid-parallel algorithm are also proposed for the types of operations typical of atmospheric inverse problems, which further reduce the cost of sparse matrix multiplication by yielding only upper triangular and/or dense matrices.

Reference Paper:

Technical Note: Improving the computational efficiency of sparse matrix multiplication in linear atmospheric inverse problems 

Yadav, V. and Michalak, A. M.: Technical Note: Improving the computational efficiency of sparse matrix multiplication in linear atmospheric inverse problems, Geosci. Model Dev. Discuss., https://doi.org/10.5194/gmd-2016-204.

Downloadable from: https://www.geosci-model-dev-discuss.net/gmd-2016-204/gmd-2016-204.pdf

Code Description:

The C++ (mixed C and C++) code demonstrates application of the operations mentioned under subheading Primary algorithm. The code is extremely efficient and uses openmp. Note filereading_main.cpp is not required to compile the code. However it was designed to test the performance of the code for large sparse matrices.
