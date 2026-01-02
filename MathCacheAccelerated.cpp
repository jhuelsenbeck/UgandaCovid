#include "MathCacheAccelerated.hpp"
#include <cmath>
#include <cstring>
#include <algorithm>

// -------------------------------------------------------------------
// Padé[13/13] coefficients for matrix exponential
// -------------------------------------------------------------------
// These are the coefficients b_k for the (13,13) Padé approximant:
//   exp(A) ≈ [p_13(A)]^{-1} * q_13(A)
// where p and q are degree-13 polynomials.
//
// Reference: Higham, "The Scaling and Squaring Method for the Matrix
// Exponential Revisited", SIAM J. Matrix Anal. Appl., 2005.
// -------------------------------------------------------------------

static const double PADE_B[14] = {
    64764752532480000.0,    // b_0
    32382376266240000.0,    // b_1
    7771770303897600.0,     // b_2
    1187353796428800.0,     // b_3
    129060195264000.0,      // b_4
    10559470521600.0,       // b_5
    670442572800.0,         // b_6
    33522128640.0,          // b_7
    1323241920.0,           // b_8
    40840800.0,             // b_9
    960960.0,               // b_10
    16380.0,                // b_11
    182.0,                  // b_12
    1.0                     // b_13
};

// Theta_13: the matrix norm threshold for Padé[13/13]
// If ||A|| > theta_13, we need to scale A down first
static const double THETA_13 = 5.371920351148152;



MathCacheAccelerated::MathCacheAccelerated(void) : MathCache() {

    workSize = 0;
}

MathCacheAccelerated::~MathCacheAccelerated(void) {

}

void MathCacheAccelerated::ensureWorkBuffers(size_t n) {

    size_t nn = n * n;
    if (workSize != nn) {
        workA.resize(nn);
        workA2.resize(nn);
        workA4.resize(nn);
        workA6.resize(nn);
        workU.resize(nn);
        workV.resize(nn);
        workTemp1.resize(nn);
        workTemp2.resize(nn);
        workTemp1T.resize(nn);
        workTemp2T.resize(nn);
        workPivot.resize(n);
        workSize = nn;
    }
}

// -------------------------------------------------------------------
// Matrix multiplication - uses BLAS when available
// -------------------------------------------------------------------

void MathCacheAccelerated::optimizedMultiply(double* A, double* B, double* C, int n) {

#if defined(USE_ACCELERATE) || defined(USE_CBLAS)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A, n, B, n, 0.0, C, n);
#else
    portableMultiply(A, B, C, n);
#endif
}

void MathCacheAccelerated::portableMultiply(double* A, double* B, double* C, int n) {

    // Cache-friendly blocked matrix multiplication
    // Block size chosen for typical L1 cache
    const int BLOCK = 32;
    
    std::fill(C, C + n*n, 0.0);
    
    for (int ii = 0; ii < n; ii += BLOCK) {
        int imax = std::min(ii + BLOCK, n);
        for (int jj = 0; jj < n; jj += BLOCK) {
            int jmax = std::min(jj + BLOCK, n);
            for (int kk = 0; kk < n; kk += BLOCK) {
                int kmax = std::min(kk + BLOCK, n);
                
                for (int i = ii; i < imax; i++) {
                    for (int k = kk; k < kmax; k++) {
                        double aik = A[i * n + k];
                        for (int j = jj; j < jmax; j++) {
                            C[i * n + j] += aik * B[k * n + j];
                        }
                    }
                }
            }
        }
    }
}

// -------------------------------------------------------------------
// Portable LU solve fallback
// -------------------------------------------------------------------

void MathCacheAccelerated::portableLUSolve(double* A, double* B, int n) {

    // In-place LU decomposition with partial pivoting
    std::vector<int> perm(n);
    for (int i = 0; i < n; i++)
        perm[i] = i;
    
    // LU factorization
    for (int k = 0; k < n - 1; k++) {
        // Find pivot
        int maxRow = k;
        double maxVal = std::fabs(A[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double val = std::fabs(A[i * n + k]);
            if (val > maxVal) {
                maxVal = val;
                maxRow = i;
            }
        }
        
        // Swap rows
        if (maxRow != k) {
            std::swap(perm[k], perm[maxRow]);
            for (int j = 0; j < n; j++) {
                std::swap(A[k * n + j], A[maxRow * n + j]);
            }
        }
        
        // Eliminate
        if (std::fabs(A[k * n + k]) > 1e-15) {
            for (int i = k + 1; i < n; i++) {
                double factor = A[i * n + k] / A[k * n + k];
                A[i * n + k] = factor;
                for (int j = k + 1; j < n; j++) {
                    A[i * n + j] -= factor * A[k * n + j];
                }
            }
        }
    }
    
    // Solve for each column of B
    for (int col = 0; col < n; col++) {
        // Apply permutation and forward substitution (Ly = Pb)
        std::vector<double> y(n);
        for (int i = 0; i < n; i++) {
            y[i] = B[perm[i] * n + col];
            for (int j = 0; j < i; j++) {
                y[i] -= A[i * n + j] * y[j];
            }
        }
        
        // Back substitution (Ux = y)
        for (int i = n - 1; i >= 0; i--) {
            double sum = y[i];
            for (int j = i + 1; j < n; j++) {
                sum -= A[i * n + j] * B[j * n + col];
            }
            B[i * n + col] = sum / A[i * n + i];
        }
    }
}

double MathCacheAccelerated::infinityNorm(double* A, size_t n) {

    double maxSum = 0.0;
    for (size_t i = 0; i < n; i++) {
        double rowSum = 0.0;
        for (size_t j = 0; j < n; j++)
            rowSum += std::fabs(A[i * n + j]);
        if (rowSum > maxSum)
            maxSum = rowSum;
    }
    return maxSum;
}

void MathCacheAccelerated::computeMatrixExponential(const DoubleMatrix& Q, double t, DoubleMatrix& P) {

    size_t n = Q.getNumRows();
    size_t nn = n * n;
    int N = static_cast<int>(n);
    
    // Ensure working buffers are allocated
    ensureWorkBuffers(n);
    
    // Use pointers for cleaner code
    double* A  = workA.data();
    double* A2 = workA2.data();
    double* A4 = workA4.data();
    double* A6 = workA6.data();
    double* U  = workU.data();
    double* V  = workV.data();
    double* temp1  = workTemp1.data();
    double* temp2  = workTemp2.data();
    double* temp1T = workTemp1T.data();
    double* temp2T = workTemp2T.data();
    
    // A = Q * t
    const double* qPtr = Q.begin();
    for (size_t i = 0; i < nn; i++)
        A[i] = qPtr[i] * t;
    
    // Compute infinity norm of A
    double normA = infinityNorm(A, n);
    
    // Scaling: if ||A|| > theta_13, scale A down by power of 2
    int s = 0;
    if (normA > THETA_13) {
        s = static_cast<int>(std::ceil(std::log2(normA / THETA_13)));
        double scale = 1.0 / (1L << s);
#if defined(USE_ACCELERATE) || defined(USE_CBLAS)
        cblas_dscal(static_cast<int>(nn), scale, A, 1);
#else
        for (size_t i = 0; i < nn; i++)
            A[i] *= scale;
#endif
    }
    
    // Compute A^2, A^4, A^6
    optimizedMultiply(A, A, A2, N);
    optimizedMultiply(A2, A2, A4, N);
    optimizedMultiply(A2, A4, A6, N);
    
    // -------------------------------------------------------------------
    // Build U = A * (A6*(b13*A6 + b11*A4 + b9*A2 + b7*I) + b5*A4 + b3*A2 + b1*I)
    // -------------------------------------------------------------------
    
    // temp1 = b7*I + b9*A2 + b11*A4 + b13*A6
    std::fill(temp1, temp1 + nn, 0.0);
    for (size_t i = 0; i < n; i++)
        temp1[i * n + i] = PADE_B[7];
    for (size_t i = 0; i < nn; i++)
        temp1[i] += PADE_B[9] * A2[i] + PADE_B[11] * A4[i] + PADE_B[13] * A6[i];
    
    // temp2 = A6 * temp1
    optimizedMultiply(A6, temp1, temp2, N);
    
    // temp2 += b5*A4 + b3*A2 + b1*I
    for (size_t i = 0; i < nn; i++)
        temp2[i] += PADE_B[5] * A4[i] + PADE_B[3] * A2[i];
    for (size_t i = 0; i < n; i++)
        temp2[i * n + i] += PADE_B[1];
    
    // U = A * temp2
    optimizedMultiply(A, temp2, U, N);
    
    // -------------------------------------------------------------------
    // Build V = A6*(b12*A6 + b10*A4 + b8*A2 + b6*I) + b4*A4 + b2*A2 + b0*I
    // -------------------------------------------------------------------
    
    // temp1 = b6*I + b8*A2 + b10*A4 + b12*A6
    std::fill(temp1, temp1 + nn, 0.0);
    for (size_t i = 0; i < n; i++)
        temp1[i * n + i] = PADE_B[6];
    for (size_t i = 0; i < nn; i++)
        temp1[i] += PADE_B[8] * A2[i] + PADE_B[10] * A4[i] + PADE_B[12] * A6[i];
    
    // V = A6 * temp1
    optimizedMultiply(A6, temp1, V, N);
    
    // V += b4*A4 + b2*A2 + b0*I
    for (size_t i = 0; i < nn; i++)
        V[i] += PADE_B[4] * A4[i] + PADE_B[2] * A2[i];
    for (size_t i = 0; i < n; i++)
        V[i * n + i] += PADE_B[0];
    
    // -------------------------------------------------------------------
    // Solve (V - U) * P = (V + U) for P
    // -------------------------------------------------------------------
    
    // temp1 = V - U (coefficient matrix)
    // temp2 = V + U (right-hand side)
    for (size_t i = 0; i < nn; i++) {
        temp1[i] = V[i] - U[i];
        temp2[i] = V[i] + U[i];
    }
    
#if defined(USE_ACCELERATE)
    // LAPACK expects column-major, so transpose both matrices
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            temp1T[j * n + i] = temp1[i * n + j];
            temp2T[j * n + i] = temp2[i * n + j];
        }
    }
    
    lapack_int info;
    lapack_int LN = static_cast<lapack_int>(n);
    
    dgetrf_(&LN, &LN, temp1T, &LN, workPivot.data(), &info);
    if (info == 0)
        dgetrs_("N", &LN, &LN, temp1T, &LN, workPivot.data(), temp2T, &LN, &info);
    
    // Transpose solution back to row-major into P
    double* pPtr = P.begin();
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            pPtr[i * n + j] = temp2T[j * n + i];
        }
    }
    
#elif defined(USE_CBLAS)
    // OpenBLAS/system LAPACK path
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            temp1T[j * n + i] = temp1[i * n + j];
            temp2T[j * n + i] = temp2[i * n + j];
        }
    }
    
    int info;
    int LN = static_cast<int>(n);
    std::vector<int> ipiv(n);
    
    dgetrf_(&LN, &LN, temp1T, &LN, ipiv.data(), &info);
    if (info == 0)
        dgetrs_("N", &LN, &LN, temp1T, &LN, ipiv.data(), temp2T, &LN, &info);
    
    double* pPtr = P.begin();
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            pPtr[i * n + j] = temp2T[j * n + i];
        }
    }
    
#else
    // Portable fallback - solve in place
    portableLUSolve(temp1, temp2, N);
    memcpy(P.begin(), temp2, nn * sizeof(double));
#endif
    
    // -------------------------------------------------------------------
    // Squaring phase: P = P^(2^s)
    // -------------------------------------------------------------------
    
    double* pOut = P.begin();
    for (int sq = 0; sq < s; sq++) {
        optimizedMultiply(pOut, pOut, temp1, N);
        memcpy(pOut, temp1, nn * sizeof(double));
    }
    
    // -------------------------------------------------------------------
    // Cleanup: ensure all entries are non-negative
    // (numerical errors can produce tiny negative values)
    // -------------------------------------------------------------------
    
    for (size_t i = 0; i < nn; i++)
        pOut[i] = std::fabs(pOut[i]);
}
