#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "main.h"

#define MAX_FILENAME_SIZE 256
#define MAX_NUM_LENGTH 100


/* This function checks the number of input parameters to the program to make 
   sure it is correct. If the number of input parameters is incorrect, it 
   prints out a message on how to properly use the program.
   input parameters:
       int    argc
       char** argv 
    return parameters:
       none
 */
void usage(int argc, char** argv)
{
    if(argc < 4) {
        fprintf(stderr, "usage: %s <matrix> <vector> <result>\n", argv[0]);
        exit(EXIT_FAILURE);
    } 
}


/* This function prints out information about a sparse matrix
   input parameters:
       char*       fileName    name of the sparse matrix file
       MM_typecode matcode     matrix information
       int         m           # of rows
       int         n           # of columns
       int         nnz         # of non-zeros
   return paramters:
       none
 */
void print_matrix_info(char* fileName, MM_typecode matcode, 
                       int m, int n, int nnz)
{
    fprintf(stdout, "-----------------------------------------------------\n");
    fprintf(stdout, "Matrix name:     %s\n", fileName);
    fprintf(stdout, "Matrix size:     %d x %d => %d\n", m, n, nnz);
    fprintf(stdout, "-----------------------------------------------------\n");
    fprintf(stdout, "Is matrix:       %d\n", mm_is_matrix(matcode));
    fprintf(stdout, "Is sparse:       %d\n", mm_is_sparse(matcode));
    fprintf(stdout, "-----------------------------------------------------\n");
    fprintf(stdout, "Is complex:      %d\n", mm_is_complex(matcode));
    fprintf(stdout, "Is real:         %d\n", mm_is_real(matcode));
    fprintf(stdout, "Is integer:      %d\n", mm_is_integer(matcode));
    fprintf(stdout, "Is pattern only: %d\n", mm_is_pattern(matcode));
    fprintf(stdout, "-----------------------------------------------------\n");
    fprintf(stdout, "Is general:      %d\n", mm_is_general(matcode));
    fprintf(stdout, "Is symmetric:    %d\n", mm_is_symmetric(matcode));
    fprintf(stdout, "Is skewed:       %d\n", mm_is_skew(matcode));
    fprintf(stdout, "Is hermitian:    %d\n", mm_is_hermitian(matcode));
    fprintf(stdout, "-----------------------------------------------------\n");

}


/* This function checks the return value from the matrix read function, 
   mm_read_mtx_crd(), and provides descriptive information.
   input parameters:
       int ret    return value from the mm_read_mtx_crd() function
   return paramters:
       none
 */
void check_mm_ret(int ret)
{
    switch(ret)
    {
        case MM_COULD_NOT_READ_FILE:
            fprintf(stderr, "Error reading file.\n");
            exit(EXIT_FAILURE);
            break;
        case MM_PREMATURE_EOF:
            fprintf(stderr, "Premature EOF (not enough values in a line).\n");
            exit(EXIT_FAILURE);
            break;
        case MM_NOT_MTX:
            fprintf(stderr, "Not Matrix Market format.\n");
            exit(EXIT_FAILURE);
            break;
        case MM_NO_HEADER:
            fprintf(stderr, "No header information.\n");
            exit(EXIT_FAILURE);
            break;
        case MM_UNSUPPORTED_TYPE:
            fprintf(stderr, "Unsupported type (not a matrix).\n");
            exit(EXIT_FAILURE);
            break;
        case MM_LINE_TOO_LONG:
            fprintf(stderr, "Too many values in a line.\n");
            exit(EXIT_FAILURE);
            break;
        case MM_COULD_NOT_WRITE_FILE:
            fprintf(stderr, "Error writing to a file.\n");
            exit(EXIT_FAILURE);
            break;
        case 0:
            fprintf(stdout, "file loaded.\n");
            break;
        default:
            fprintf(stdout, "Error - should not be here.\n");
            exit(EXIT_FAILURE);
            break;

    }
}


/* This function reads information about a sparse matrix using the 
   mm_read_banner() function and prints out information using the
   print_matrix_info() function.
   input parameters:
       char*       fileName    name of the sparse matrix file
   return paramters:
       none
 */
void read_info(char* fileName)
{
    FILE* fp;
    MM_typecode matcode;
    int m;
    int n;
    int nnz;

    if((fp = fopen(fileName, "r")) == NULL) {
        fprintf(stderr, "Error opening file: %s\n", fileName);
        exit(EXIT_FAILURE);
    }

    if(mm_read_banner(fp, &matcode) != 0)
    {
        fprintf(stderr, "Error processing Matrix Market banner.\n");
        exit(EXIT_FAILURE);
    } 

    if(mm_read_mtx_crd_size(fp, &m, &n, &nnz) != 0) {
        fprintf(stderr, "Error reading size.\n");
        exit(EXIT_FAILURE);
    }

    print_matrix_info(fileName, matcode, m, n, nnz);

    fclose(fp);
}


/* This function coverts a sparse matrix stored in COO format to CSR.
   input parameters:
       these are 'consumed' by this function
       int*           row_ind 		row index for the non-zeros in COO
       int*           col_ind		column index for the non-zeros in COO
       double*        val		values for the non-zeros in COO
       int            m			# of rows in the matrix
       int            n			# of columns in the matrix
       int            nnz		# of non-zeros in the matrix

       these are 'produced' by this function
       unsigned int** csr_row_ptr	row pointers to csr_col_ind and 
                                        csr_vals in CSR 
       unsigned int** csr_col_ind	column index for the non-zeros in CSR
       double**       csr_vals		values for the non-zeros in CSR
   return parameters:
       none
 */
void convert_coo_to_csr(int* row_ind, int* col_ind, double* val, 
                        int m, int n, int nnz,
                        unsigned int** csr_row_ptr, unsigned int** csr_col_ind,
                        double** csr_vals)
{
    // matrix row starts at 1 so m+1 (not sure why we +1)
    // maybe it's to fill at the exact same index as the coords
    // ignoring the whole 0 row because there's no 0 in the COO
    unsigned int *trow = (unsigned int *)malloc(sizeof(unsigned int) * (m+1));
    unsigned int *tcol = (unsigned int *)malloc(sizeof(unsigned int) * nnz);
    double *tval = (double *)malloc(sizeof(double) * nnz);
    // histogram
    int *buckets = (int *)malloc(sizeof(int) * m);

    if (buckets == NULL || trow == NULL || tcol == NULL || tval == NULL) {
        printf("problem with mem allocation");
        exit(EXIT_FAILURE);
    }
    memset(buckets, 0, sizeof(int));
    for (int i = 0; i < nnz; i++) {
        // there is no row 0 so buckets[0] will always equal 0
        buckets[row_ind[i] - 1] += 1;
    }
    // prefix sum
    for (int i = 1; i < m; i++) {
        buckets[i] = buckets[i] + buckets[i - 1];
    }
    // puts the prefix sum into trow with +1 so it has the 
    // same row number and index
    trow[0] = 0;
    for (int i = 0; i < m; i++) {
        // buckets[0] == 0, so now trow[0], trow[1] = 0, trow[2] = 13 ...
        trow[i+1] = (unsigned int)buckets[i];
    }
    /*
    trow has all the rows start endpointers in it
    loop over nnz and place index in pointed to spot
    offset of each row i kept in tmp_row
    we add 1 to the current row pointer to get the offset for the next
    time that row is found. We need the initial pointers to get the initial
    location to start counting each row from.
    */
    unsigned int *tmp_row = (unsigned int *)malloc(sizeof(unsigned int) * m);
    memcpy(tmp_row, trow, sizeof(unsigned int) * m);
    for (int i = 0; i < nnz; i++) {
        // col_ind[i] -1 bc col also starts at 1 in matrix
        // tmp_row/trow has 0 at index 0, which corresponds to row1
        tcol[tmp_row[row_ind[i] - 1]] = col_ind[i] - 1;
        tval[tmp_row[row_ind[i] - 1]] = val[i];
        // cumulating offset
        tmp_row[row_ind[i]-1]++;
    }

    *csr_row_ptr = trow;
    *csr_col_ind = tcol;
    *csr_vals = tval;
    free(tmp_row);
    free(buckets);
}
// set args test2/A.mtx test2/x.mtx test2/ans.mtx
// set args test2/A.mtx test2/x.mtx test2/liams.mtx


/* This function reads in a vector from a text file, similar in format to
   the Matrix Market format.
   The first line contains the number of elements in the vector.
   The rest of the file contains the values in the vector, one element per row.
   input parameters:
       char*    fileName	Name of the file containing the vector
       double** vector		Array that will contain the vector
       int*     vecSize		Integer variable that will contain the size of
                                the vector
   return parameters:
       none
 */
// vector and vecsize are initialized within this function
void read_vector(char* fileName, double** vector, int* vecSize)
{
    // reads vecSize for mallocing vector
    FILE* fp = fopen(fileName, "r");
    assert(fp);
    char line[MAX_NUM_LENGTH]; // stores number of rows in vector
    fgets(line, MAX_NUM_LENGTH, fp);
    fclose(fp);

    unsigned int vector_size = atoi(line);
    double* vector_ = (double*) malloc(sizeof(double) * vector_size);

    fp = fopen(fileName, "r");
    assert(fp); 
    // first read the first line again to move file pointer to first row
    fgets(line, MAX_NUM_LENGTH, fp);

    unsigned int index = 0;
    while(fgets(line, MAX_NUM_LENGTH, fp) != NULL) {
        // atof: string -> floating point
        vector_[index] = atof(line); 
        index++;
    }

    fclose(fp);
    assert(index == vector_size);

    *vector = vector_;
    *vecSize = vector_size;
}

/* This function calculates the sparse matrix-vector multiply from the matrix
   in CSR format (i.e., csr_row_ptr, csr_col_ind, and csr_vals) and the vector 
   in an array (i.e., vector_x). It stores the result in another array (i.e.,
   res)
   input parameters:
       these are 'consumed' by this function
       unsigned int** csr_row_ptr	row pointers to csr_col_ind and 
                                        csr_vals in CSR 
       unsigned int** csr_col_ind	column index for the non-zeros in CSR
       double**       csr_vals		values for the non-zeros in CSR
       int            m			# of rows in the matrix
       int            n			# of columns in the matrix
       int            nnz		# of non-zeros in the matrix
       double         vector_x		input vector

       these are 'produced' by this function
       double*        res		Result of SpMV. res = A * x, where
                                        A is stored in CSR format and x is 
                                        stored in vector_x
   return parameters:
       none

 */
// multiplying a csr, A, by a vector, vector_x
void spmv(unsigned int* csr_row_ptr, unsigned int* csr_col_ind, 
          double* csr_vals, int m, int n, int nnz, 
          double* vector_x, double* res)
{
    // loop over rows
    for (int i = 0; i < m; i++) {
        // loop over cols
        for (int j = csr_row_ptr[i]; j < csr_row_ptr[i+1]; j++) {
            res[i] += (double)(csr_vals[j]) * (double)(vector_x[csr_col_ind[j]]);
        }
    }
}


/* This function stores a vector in a text file, similar in format to
   the Matrix Market format.
   The first line contains the number of elements in the vector.
   The rest of the file contains the values in the vector, one element per row.
   input parameters:
       char*    fileName	Name of the file that will contain the vector
       double** res 		Array that contains the vector
       int*     m		Integer variable that contains the size of
                                the vector
   return parameters:
       none
 */
void store_result(char *fileName, double* res, int m)
{
    FILE* fp = fopen(fileName, "w");
    assert(fp);

    fprintf(fp, "%d\n", m);
    for(int i = 0; i < m; i++) {
        fprintf(fp, "%0.10f\n", res[i]);
    }

    fclose(fp);
}

/* This program first reads in a sparse matrix stored in matrix market format 
   (mtx). It generates a set of arrays - row_ind, col_ind, and val, which stores
   the row/column index and the value for the non-zero elements in the matrix, 
   respectively. This is also known as the co-ordinate format.

   Then, it should convert this matrix stored in co-ordinate format to the
   compressed sparse row (CSR) format.

   Then, finally it should use the CSR format to multiply the matrix with a 
   vector (i.e., calculate the sparse matrix vector multiply, or SpMV).

   The resulting vector should then be stored in a file, one value per line, 
   whose name was specified as an input to the program.
 */
int main(int argc, char** argv)
{
    usage(argc, argv);

    // Read the sparse matrix file name
    char matrixName[MAX_FILENAME_SIZE];
    strcpy(matrixName, argv[1]);
    read_info(matrixName);

    // Read the sparse matrix and store it in row_ind, col_ind, and val,
    // also known as co-ordinate format.
    int ret;
    MM_typecode matcode;
    int m;
    int n;
    int nnz;
    int *row_ind; // pointer to csr row
    int *col_ind; // pointer to csr col
    double *val; // points to csr val

    fprintf(stdout, "Matrix file name: %s ... ", matrixName);
    /* 
     mm_read_mtx_crd is a fucntion provided in mmio.c that reads in a sparse
     matrix in Matrix Market format and stores the matrix in COO format.
     m - # of rows
     n - # of columns
     nnz - number of non-zeroes
     row_ind - array of row indices for the non-zeros
     col_ind - array of column indices for the non-zeros
     val     - array of values for the non-nzeros
     matcode - return value for the function

     The first non-zero's row and column indices are stored in row_ind[0], and
     col_ind[0], respectively, and the value of the non-zero is stored in
     val[0].
     Therefore, the size of these arrays are equal to nnz.
     */
    // returns SpMV in COO format
    ret = mm_read_mtx_crd(matrixName, &m, &n, &nnz, &row_ind, &col_ind, &val, 
                          &matcode);
    // prints errors if any exist
    check_mm_ret(ret);

    // Convert co-ordinate format to CSR format
    fprintf(stdout, "Converting COO to CSR...");
    // pointers to start of csr row and col respectively
    // rememebr csr_row_ptr holds index of csr col where new row starts
    unsigned int* csr_row_ptr = NULL; 
    unsigned int* csr_col_ind = NULL;  
    // pointer to csr val
    double* csr_vals = NULL; 
    convert_coo_to_csr(row_ind, col_ind, val, m, n, nnz,
                       &csr_row_ptr, &csr_col_ind, &csr_vals);
    fprintf(stdout, "done\n");

    // Load the vector file
    char vectorName[MAX_FILENAME_SIZE];
    strcpy(vectorName, argv[2]);
    fprintf(stdout, "Vector file name: %s ... ", vectorName);
    double* vector_x;
    unsigned int vector_size;
    read_vector(vectorName, &vector_x, &vector_size);
    assert(n == vector_size);
    fprintf(stdout, "file loaded\n");

    // Calculate SpMV
    double *res = (double*) malloc(sizeof(double) * m);;
    assert(res);
    fprintf(stdout, "Calculating SpMV ... ");
    spmv(csr_row_ptr, csr_col_ind, csr_vals, m, n, nnz, vector_x, res);
    fprintf(stdout, "done\n");

    // Store the calculated vector in a file, one element per line.
    char resName[MAX_FILENAME_SIZE];
    strcpy(resName, argv[3]); 
    fprintf(stdout, "Result file name: %s ... ", resName);
    store_result(resName, res, m);
    fprintf(stdout, "file saved\n");

    // Free memory
    free(csr_row_ptr);
    free(csr_col_ind);
    free(csr_vals);

    free(vector_x);
    free(res);

    free(row_ind);
    free(col_ind);
    free(val);

    return 0;
}
