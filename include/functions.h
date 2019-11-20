
#ifndef FUNCTIONS_INCLUDED
#define FUNCTIONS_INCLUDED

void destroy(struct sparsemat *matrix);

void destroy(struct iarray *matrix);

void destroy(struct darray *matrix);

void modifyalloc(struct sparsemat *matrix, int imemSize);

struct sparsemat *create(int rows, int cols, int nzmax);

struct darray *create(int rows, int cols);

struct iarray limits(int tcov_rows, int numprocs);

void sparsework_nosym(const struct sparsemat * const matrixa, const struct sparsemat * const matrixb, struct sparsemat * const matrixc, const int startIndex, const int endIndex, int memIncrease);

void sparsework_sym(const struct sparsemat * const matrixa, const struct sparsemat * const matrixb, struct sparsemat * const matrixc, const int startIndex, const int endIndex, int memIncrease);

void dense_sym(const struct sparsemat * const matrixa, const struct sparsemat * const matrixb, struct darray * const matrixc);

void dense_nosym(const struct sparsemat * const matrixa, const struct sparsemat * const matrixb, struct darray * const matrixc);

void sparse_nosym(const struct sparsemat * const  matrixA, const struct sparsemat* const  matrixB, struct sparsemat *  const matrixC, int imemSize);

void sparse_sym(const struct sparsemat * const  matrixA, const struct sparsemat* const  matrixB, struct sparsemat *  const matrixC, int imemSize);

void print_array_1d(double * __restrict Array, int length);

void printMatForm(double * __restrict Array, const int rows, const int cols, int width, const int precision);

void print_array_1d(int * __restrict Array, int length);

void printMatForm(int * __restrict Array, const int rows, const int cols, int width, const int precision);

void sparsemult_single(const struct sparsemat * const matrixa, const struct sparsemat * const matrixb, struct sparsemat * matrixc, int imemsize);

#endif // Memory functions included
