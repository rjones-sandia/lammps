#ifndef VOIGT_OPERATIONS_H
#define VOIGT_OPERATIONS_H
#include "MatrixDef.h"
#include "MatrixLibrary.h"

// LAMMPS stress/atom order 11,22,33,12,13,23 and Output
// Voigt                    11,22,33,23,13,12

// Voigt indexing puts a symmetric 3x3 matrix into a 
// vector form: [11 22 33 12 13 23]
//
// matrix form: [[ 0 5 4 ] 
//               [ 5 1 3 ]
//               [ 4 3 2 ]]
//
// unsymmetric Voigt
// vector form: [11 22 33 12 13 23 21 31 32]
//
// matrix form: [[ 0 3 4 ] 
//               [ 8 1 4 ]
//               [ 7 6 2 ]]
  
// Matrix order for Output
// vector form: [11 12 13 21 22 23 31 32 33]
//
// matrix form: [[ 0 1 2 ] 
//               [ 3 4 5 ]
//               [ 6 7 8 ]]

namespace voigt3 {
  const int voigt_idx1[] = {0,1,2,0,0,1,1,2,2}; // first  packed voigt index
  const int voigt_idx2[] = {0,1,2,1,2,2,0,0,1}; // second packed voigt index
  //* Returns the trace of a 3x3 matrix in symmetric voigt form.
  inline double tr(const DENS_VEC &A)
  {
      return A(0) + A(1) + A(2);
  }
  //* Returns the identify matrix in voigt form, optionally scaled by a factor.
  inline DENS_VEC eye(double scale=1.0)
  {
    DENS_VEC I(6);
    I(0) = I(1) = I(2) = scale;
    I(3) = I(4) = I(5) = 0.;
    return I;
  } 
  //* Computes a symmetric matrix-matrix product
  //* Inputs 6-length vectors A, B
  inline DENS_VEC dsymm(const DENS_VEC &A, const DENS_VEC &B)
  {
    DENS_VEC C(6,false); 
    C(0) = A(0)*B(0)+A(3)*B(3)+A(4)*B(4); // 11*11+12*21+13*31
    C(1) = A(3)*B(3)+A(1)*B(1)+A(5)*B(5); // 21*12+22*22+23*32
    C(2) = A(4)*B(4)+A(5)*B(5)+A(2)*B(2); // 31*13+32*23+33*33
    C(3) = A(0)*B(3)+A(3)*B(1)+A(4)*B(5); // 11*12+12*22+13*32
    C(4) = A(0)*B(4)+A(3)*B(5)+A(4)*B(2); // 11*13+12*23+13*33
    C(5) = A(3)*B(4)+A(1)*B(5)+A(5)*B(2); // 21*13+22*23+23*33
    return C;
  }
  //* Computes the determinant of a 3x3 matrix in symmetric voigt form.
  inline double det(const DENS_VEC &A)
  {
    return A(0) * (A(1)*A(2)-A(5)*A(5))  // 11*(22*33-23*32)
          -A(3) * (A(3)*A(2)-A(5)*A(4))  // 12*(21*33-31*32)
          +A(4) * (A(3)*A(5)-A(1)*A(4)); // 13*(21*32-22*31)
  }
  //* Returns the derivative of C*C in voigt notation.
  inline DENS_MAT derivative_of_square(const DENS_VEC &C)
  {
    DENS_MAT D(6,6);
    D(0,0)=2.0*C(0);  D(0,1)=0.0;       D(0,2)=0.0; 
    D(1,0)=0.0;       D(1,1)=2.0*C(1);  D(1,2)=0.0; 
    D(2,0)=0.0;       D(2,1)=0.0;       D(2,2)=2.0*C(2); 

    D(0,3)=0.0;       D(0,4)=2.0*C(4);  D(0,5)=2.0*C(5);
    D(1,3)=2.0*C(3);  D(1,4)=0.0;       D(1,5)=2.0*C(5);
    D(2,3)=2.0*C(3);  D(2,4)=2.0*C(4);  D(2,5)=0.0;

    D(3,0)=0.0;       D(3,1)=C(3);      D(3,2)=C(3); 
    D(4,0)=C(4);      D(4,1)=0.0;       D(4,2)=C(4); 
    D(5,0)=C(5);      D(5,1)=C(5);      D(5,2)=0.0; 
   
    D(3,3)=C(0)+C(1); D(3,4)=C(5);      D(3,5)=C(4);
    D(4,3)=C(5);      D(4,4)=C(0)+C(2); D(4,5)=C(3);
    D(5,3)=C(4);      D(5,4)=C(3);      D(5,5)=C(1)+C(2); 
    return D;
  }
  //* Computes the inverse of a 3x3 matrix in symmetric voigt form.
  inline DENS_VEC inv(const DENS_VEC &A)
  {
    DENS_VEC B(6,false);
    const double inv_det = 1.0/det(A);
    B(0) = (A(1)*A(2)-A(5)*A(5))*inv_det; // 22*33-23*32
    B(1) = (A(0)*A(2)-A(4)*A(4))*inv_det; // 11*33-13*31
    B(2) = (A(0)*A(1)-A(3)*A(3))*inv_det; // 11*22-12*32
    B(3) = (A(4)*A(5)-A(3)*A(2))*inv_det; // 31*23-21*33
    B(4) = (A(3)*A(5)-A(4)*A(1))*inv_det; // 21*32-31*22
    B(5) = (A(4)*A(3)-A(0)*A(5))*inv_det; // 31*12-11*32
    return B;
  }
  //* Returns the voigt form of a symmetric matrix.
  // consistent with voigt_idx1,2
  inline DENS_VEC to_voigt(const DENS_MAT &C)
  {
    DENS_VEC B(6,false);
    B(0)=C(0,0);
    B(1)=C(1,1);
    B(2)=C(2,2);
    B(3)=C(0,1); // take upper triangle entries
    B(4)=C(0,2); 
    B(5)=C(1,2);  
    return B; 
  }
  //* Returns a symmetric matrix form a voigt form.
  // consistent with voigt_idx1,2
  inline DENS_MAT from_voigt(const DENS_VEC &B)
  {
    DENS_MAT C(3,3,false);
    C(0,0)=B(0);  C(0,1)=B(3);  C(0,2)=B(4);
    C(1,0)=B(3);  C(1,1)=B(1);  C(1,2)=B(5);
    C(2,0)=B(4);  C(2,1)=B(5);  C(2,2)=B(2);
    return C;
  }
  //* Returns the voigt form of an unsymmetric matrix.
  // consistent with voigt_idx1,2
  inline DENS_VEC to_voigt_unsymmetric(const DENS_MAT &C)
  {
    DENS_VEC B(9,false);
    B(0)=C(0,0);
    B(1)=C(1,1);
    B(2)=C(2,2);
    B(3)=C(0,1); // upper triangle entries
    B(4)=C(0,2); 
    B(5)=C(1,2);  
    B(6)=C(1,0); // lower triangle entries
    B(7)=C(2,0); 
    B(8)=C(2,1);  
    return B; 
  }
  inline DENS_VEC to_voigt_unsymmetric(const DENS_VEC &C)
  {
    DENS_VEC B(9,false);
    B(0)=C(0);
    B(1)=C(1);
    B(2)=C(2);
    B(3)=C(3); // upper triangle entries
    B(4)=C(4); 
    B(5)=C(5);  
    B(6)=C(3); // lower triangle entries
    B(7)=C(4); 
    B(8)=C(5);  
    return B; 
  }
  
  //* Returns a symmetric matrix form a voigt form.
  // consistent with voigt_idx1,2
  inline DENS_MAT from_voigt_unsymmetric(const DENS_VEC &B)
  {
    DENS_MAT C(3,3,false);
    C(0,0)=B(0);  C(0,1)=B(3);  C(0,2)=B(4);
    C(1,0)=B(6);  C(1,1)=B(1);  C(1,2)=B(5);
    C(2,0)=B(7);  C(2,1)=B(8);  C(2,2)=B(2);
    return C;
  }

  //* Converts vector to DENS_MAT_VEC 
  inline void vector_to_dens_mat_vec(const DENS_MAT & IN, DENS_MAT_VEC & OUT)
  {
    for (int i=0; i<IN.nRows(); i++) {
      for (int j=0; j<3; j++) {
        for (DENS_MAT_VEC::size_type k=0; k<3; k++) {
          OUT[k](i,j) = IN(i,3*j+k);
         }
      }
    }
    return;
  }
  //* Converts DENS_MAT_VEC to vector
  inline void symm_dens_mat_vec_to_vector(const DENS_MAT_VEC & IN, DENS_MAT & OUT)
  {
    for (int i=0; i<IN.front().nRows(); i++) {
      for (int v=0; v<6; v++) {
        OUT(i,v) = IN[voigt_idx1[v]](i,voigt_idx2[v]);
      }
    }
    return;
  }
  inline void global_to_vector(const int i, const DENS_MAT & IN, DENS_VEC & OUT)
  {
    int ncols = IN.nCols();
    OUT.reset(ncols);
    for (int j = 0; j < ncols; ++j){ OUT(j)=IN(i,j); }
  }
  inline void global_to_matrix(const int i, const DENS_MAT & IN, DENS_MAT & OUT)
  {
    int ncols = IN.nCols();
    OUT.reset(3,3);
    if (ncols == 6) {
      OUT(0,0)=IN(i,0); OUT(0,1)=IN(i,3); OUT(0,2)=IN(i,4);  
      OUT(1,0)=IN(i,3); OUT(1,1)=IN(i,1); OUT(1,2)=IN(i,5);  
      OUT(2,0)=IN(i,4); OUT(2,1)=IN(i,5); OUT(2,2)=IN(i,2);  
    }
    else {
      OUT(0,0)=IN(i,0); OUT(0,1)=IN(i,1); OUT(0,2)=IN(i,2);  
      OUT(1,0)=IN(i,3); OUT(1,1)=IN(i,4); OUT(1,2)=IN(i,5);  
      OUT(2,0)=IN(i,6); OUT(2,1)=IN(i,7); OUT(2,2)=IN(i,8);  
    }
  }
  inline void matrix_to_global(const DENS_MAT & IN, int i, DENS_MAT & OUT)
  {
    int ncols = IN.nCols();
    if (ncols == 6) {
      OUT(i,0) = IN(0,0);
      OUT(i,1) = IN(1,1);
      OUT(i,2) = IN(2,2);
      OUT(i,3) = IN(0,1);
      OUT(i,4) = IN(0,2);
      OUT(i,5) = IN(1,2);
    }
    else {
      OUT(i,0) = IN(0,0);
      OUT(i,1) = IN(0,1);
      OUT(i,2) = IN(0,2);
      OUT(i,3) = IN(1,0);
      OUT(i,4) = IN(1,1);
      OUT(i,5) = IN(1,2);
      OUT(i,6) = IN(2,0);
      OUT(i,7) = IN(2,1);
      OUT(i,8) = IN(2,2);
    }
  }
}
//====================================================================
namespace matrix3 {
//====================================================================
  inline void add_identity(DENS_MAT & IN)
  {
    IN(0,0) += 1;
    IN(1,1) += 1;
    IN(2,2) += 1;
  }

}
#endif

