#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_zeroMat(double *Mat,int row,int col)
{
  int i,j;
  for (i = 0; i < row; i++)
    for (j = 0; j < col; j++)
        *(Mat + i*col+j) = 0;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_multVec2SkewMat(double *ptr_b, int row, double *ptr_A) //Multiply matrix to vector
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{

  if (row != 3)
    printf("Error vector should have 3 elements \n");

  int i,j;
  double A[row][row], b[row];

  ram_zeroMat(&A[0][0],row,row);

  for (i=0;i<row;i++)
		b[i] = *(ptr_b + i);

  A[0][1] = -b[2];
  A[1][0] = b[2];
  A[0][2] = b[1];
  A[2][0] = -b[1];
  A[1][2] = -b[0];
  A[2][1] = b[0];

	for (i=0;i<row;i++)
		for (j=0;j<row;j++)
			*(ptr_A + i*row+j) = A[i][j];


}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_regularizeMat(double *Mat,int row,int col, double eps)
{
  if (row != col)
  {
     printf("ERROR: Inverse currently works only on square matrices \n"); exit(1);
  }

  int i;
  for (i = 0; i < row; i++)
        *(Mat + i*col+i) += eps;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_randMat(double *Mat,int row,int col)
{
  int i,j;
  double range = 1.0 * RAND_MAX;
  for (i = 0; i < row; i++)
    for (j = 0; j < col; j++)
        *(Mat + i*col+j) = rand() / range;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_randVec(double *Vec,int row)
{
  int i;
  double range = 1.0 * RAND_MAX;
  for (i = 0; i < row; i++)
        *(Vec + i) = rand() / range;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_multMatVec(double *ptr_A, int row_A, int col_A, double *ptr_b, int row_b, double *ptr_x) //Multiply matrix to vector
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
	int i, j, p;
	double A[row_A][col_A], b[row_b], sum;

	if (col_A != row_b)
	{ printf("ERROR in multiplyMAT2VEC: Columns of A should equal rows of b \n");}// exit(1); }

	for (i=0;i<row_A;i++)
		for (j=0;j<col_A;j++)
			A[i][j] = *(ptr_A + i*col_A+j);

	for (i=0;i<row_b;i++)
		b[i] = *(ptr_b + i);

	for (i=0;i<row_A;i++)
	{
		sum = 0.0;
		for (p=0;p<col_A;p++)
			sum = sum + A[i][p]*b[p];

		*(ptr_x + i) = sum;
    }

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_multMatMat(double *ptr_A, int row_A, int col_A, double *ptr_B, int row_B, int col_B,  double *ptr_Product) //Multiply two matrices
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
	int i, j, p;
	double A[row_A][col_A], B[row_B][col_B], sum;

	if (col_A != row_B)
	{ printf("ERROR: Columns of A should equal rows of B \n"); //exit(1);
	}

	for (i=0;i<row_A;i++)
		for (j=0;j<col_A;j++)
			A[i][j] = *(ptr_A + i*col_A+j);

	for (i=0;i<row_B;i++)
		for (j=0;j<col_B;j++)
			B[i][j] = *(ptr_B + i*col_B+j);

	for (i=0;i<row_A;i++)
		for (j=0;j<col_B;j++)
		{
			sum = 0;
			for (p=0;p<col_A;p++)
				sum = sum + A[i][p]*B[p][j];

			*(ptr_Product + i*col_B+j) = sum;
		}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_initMat(double *Mat,int row,int col)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
  //initialize Matrix to 0
  int i,j;
  for (i=0;i<row;i++)
    for (j=0;j<col;j++)
      *(Mat + i*col+j) = 0; //advance one column at a time
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_mulScalarMat(double scalar,double *ptr_A, int row, int col) //multiple matrix by a scalar
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
  int i, j;
  double A[row][col];

	for (i=0;i<row;i++)
		for (j=0;j<col;j++)
		  A[i][j] = *(ptr_A + i*col+j);

	for (i=0;i<row;i++)
		for (j=0;j<col;j++)
			*(ptr_A + i*col+j) = scalar*A[i][j];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_initVec(double *Vec,int row)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
  //initialize Vector to 0
  int i;
  for (i=0;i<row;i++)
      *(Vec + i) = 0; //advance one column at a time
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_identity(double *Mat,int row) //generate identity matrix
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
	int i,j;
	for (i=0;i<row;i++)
		for (j=0;j<row;j++)
		{
			if (i==j)
				*(Mat + i*row+j) = 1.0;
			else
				*(Mat + i*row+j) = 0.0;
		}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_transpose(double *ptr_A, int row, int col,  double *ptr_At) //Take transpose of a matrix
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
	int i,j;
	double At[col][row]; //Is the transpose

	//if (row != col)
	//{ printf("ERROR: Transpose currently works only on square matrices \n"); exit(1); }

	for (i=0;i<row;i++)
		for (j=0;j<col;j++)
        {
			At[j][i] = *(ptr_A + i*col+j); //get the transpose

        }

//    for (i=0;i<col;i++)
//    {
//		for (j=0;j<row;j++)
//        {
//            printf("%f ",At[i][j]);
//        }
//        printf("\n");
//    }

	for (i=0;i<col;i++)
		for (j=0;j<row;j++)
			*(ptr_At + i*row+j) = At[i][j]; //stash back and send
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_inverse(double *ptr_A,int row, int col, double *ptr_Ainv) //invert matrix
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
  if (row != col)
  {
     printf("ERROR: Inverse currently works only on square matrices \n"); exit(1);
  }

  int i, j;
  gsl_matrix * gsl_A = gsl_matrix_alloc(row,col);
  gsl_matrix * gsl_invA = gsl_matrix_alloc(col,row);
  for (i = 0; i < row; i++)
  for (j = 0; j < col; j++)
    {
    gsl_matrix_set (gsl_A, i, j, *(ptr_A + i*col+j));
    }

    gsl_permutation *p = gsl_permutation_alloc(row);
    int s;

    // Compute the LU decomposition of this matrix
   gsl_linalg_LU_decomp(gsl_A, p, &s);
   // Compute the  inverse of the LU decomposition
   gsl_linalg_LU_invert(gsl_A, p, gsl_invA);

    for (i = 0; i < col; i++)
    for (j = 0; j < row; j++)
    {
      *(ptr_Ainv + i*row+j) = gsl_matrix_get(gsl_invA, i, j);
    }

     gsl_matrix_free (gsl_A);
     gsl_matrix_free (gsl_invA);

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_pinverse(double *ptr_A,int n_fil, int n_col, double *ptr_Ainv) //invert matrix
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
  unsigned i = 0;
	unsigned j = 0;

  double mat[n_fil][n_col];

	for (i=0;i<n_fil;i++)
		for (j=0;j<n_col;j++)
		  mat[i][j] = *(ptr_A + i*n_col+j);

	gsl_matrix * gA = gsl_matrix_alloc (n_fil, n_col);
	for (i = 0; i < n_fil; i++)
		for (j = 0; j < n_col; j++)
	   		gsl_matrix_set (gA, i, j, mat[i][j]);


	gsl_matrix * gA_t = gsl_matrix_alloc (n_col, n_fil);
	gsl_matrix_transpose_memcpy (gA_t, gA);					// Computing the transpose of gA

	gsl_matrix * U = gsl_matrix_alloc (n_col, n_fil);
	gsl_matrix * V= gsl_matrix_alloc (n_fil, n_fil);
	gsl_vector * S = gsl_vector_alloc (n_fil);


	// Computing the SVD of the transpose of A
	// The matrix 'gA_t' will contain 'U' after the function is called
	gsl_vector * work = gsl_vector_alloc (n_fil);
	gsl_linalg_SV_decomp (gA_t, V, S, work);
	gsl_vector_free(work);

	gsl_matrix_memcpy (U, gA_t);


	//Inverting S//
	//----------------------------------------------------------
	// Matrix 'S' is diagonal, so it is contained in a vector.
	// We operate to convert the vector 'S' into the matrix 'Sp'.
	//Then we invert 'Sp' to 'Spu'
	//----------------------------------------------------------
	gsl_matrix * Sp = gsl_matrix_alloc (n_fil, n_fil);
	gsl_matrix_set_zero (Sp);
	for (i = 0; i < n_fil; i++)
		gsl_matrix_set (Sp, i, i, gsl_vector_get(S, i));	// Vector 'S' to matrix 'Sp'

	gsl_permutation * p = gsl_permutation_alloc (n_fil);
	int signum;
	gsl_linalg_LU_decomp (Sp, p, &signum);				// Computing the LU decomposition

	// Compute the inverse like in the MATLAB script

	gsl_matrix * SI = gsl_matrix_calloc (n_fil, n_fil);

	for (i = 0; i < n_fil; i++) {
	  //std::cout << "S [" << i << "] = " << gsl_vector_get (S, i) << std::endl;

	  if (gsl_vector_get (S, i) > 0.0000000001)
	    gsl_matrix_set (SI, i, i, 1.0 / gsl_vector_get (S, i));
	}

	gsl_matrix * VT = gsl_matrix_alloc (n_fil, n_fil);
	gsl_matrix_transpose_memcpy (VT, V);					// Tranpose of V


	//THE PSEUDOINVERSE//
	//----------------------------------------------------------
	//Computation of the pseudoinverse of trans(A) as pinv(A) = U·inv(S).trans(V)   with trans(A) = U.S.trans(V)
	//----------------------------------------------------------
	gsl_matrix * SIpVT = gsl_matrix_alloc (n_fil, n_fil);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,				// Calculating  inv(S).trans(V)
                	1.0, SI, VT,
                	0.0, SIpVT);


	gsl_matrix * pinv = gsl_matrix_alloc (n_col, n_fil);	// Calculating  U·inv(S).trans(V)
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                	1.0, U, SIpVT,
                	0.0, pinv);

	//end THE PSEUDOINVERSE//
  for (i = 0; i < n_col; i++)
      for (j = 0; j < n_fil; j++)
        *(ptr_Ainv + i*n_fil+j) = gsl_matrix_get (pinv, i, j);
//        printf ("m(%d,%d) = %g\n", i, j,


  gsl_matrix_free(VT);
	gsl_matrix_free(SI);
	gsl_matrix_free(SIpVT);
	gsl_matrix_free(gA_t);
	gsl_matrix_free(U);
	gsl_matrix_free(gA);
	gsl_matrix_free(V);
	gsl_vector_free(S);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_display(double *Mat, int row, int col, char *name) //display the given matrix
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
    //display Matrix in MATLAB compatible format
    int i, j;

   for(i=0; i<strlen(name); i++)
   {
       printf("%c",name[i]);
   }
    printf("=[ ... \n");

	for (i=0;i<row;i++)
    {
		for (j=0;j<col;j++)
			// printf("%g ",*(Mat + i*col+j)); //advance one column at a time
      printf("%1.3f ",*(Mat + i*col+j)); //advance one column at a time
        printf("; "); printf("\n");
    }
    printf("];\n \n");
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_subtract(double *ptr_A, int row, int col, double *ptr_B, double *ptr_Sum) //add two matrices
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
    int i, j;
    double A[row][col], B[row][col];

	for (i=0;i<row;i++)
		for (j=0;j<col;j++)
		{ A[i][j] = *(ptr_A + i*col+j); B[i][j] = *(ptr_B + i*col+j);}

	for (i=0;i<row;i++)
		for (j=0;j<col;j++)
			*(ptr_Sum + i*col+j) = A[i][j] - B[i][j];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_copyMat(double *ptr_A, int row, int col, double *ptr_B) //copy one matrix into another one
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
  int i, j;

	for (i=0;i<row;i++)
		for (j=0;j<col;j++)
		   *(ptr_B + i*col+j) = *(ptr_A + i*col+j);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void ram_add(double *ptr_A, int row, int col, double *ptr_B, double *ptr_Sum) //add two matrices
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
{
    int i, j;
    double A[row][col], B[row][col];

	for (i=0;i<row;i++)
		for (j=0;j<col;j++)
		{ A[i][j] = *(ptr_A + i*col+j); B[i][j] = *(ptr_B + i*col+j);}

	for (i=0;i<row;i++)
		for (j=0;j<col;j++)
			*(ptr_Sum + i*col+j) = A[i][j] + B[i][j];
}
