#include "MathCache.hpp"
#include "Msg.hpp"

MathCache::MathCache(void) {

    matrixCount = 0; 
    arrayCount = 0;
}

MathCache::~MathCache(void) {

}

DoubleMatrix* MathCache::pushMatrix(size_t rows, size_t columns) {

    assert(matrixCount < bufferSize -1);
    auto m = matrixBuffer + matrixCount++;
    m->create(rows, columns);
    return m;
}

DoubleMatrix* MathCache::pushMatrix(size_t size) {

    assert(matrixCount < bufferSize - 1);
    auto m = matrixBuffer + matrixCount++;
    m->create(size);
    return m;
}

void MathCache::popMatrix(int n) {

    assert(matrixCount >= n);
    matrixCount -= n;
}

DoubleArray* MathCache::pushArray(size_t size) {

    assert(arrayCount < bufferSize - 1);
    auto a = arrayBuffer + arrayCount++;
    a->create(size);
    return a;
}

void MathCache::popArray() {

    assert(arrayCount > 0);
    --arrayCount;
}

void MathCache::backSubstitutionRow(DoubleMatrix& U, double* b) {

    int n = (int)U.getNumRows();
    int np1 = n + 1;

    auto bp = b + n - 1;

    auto urow  = U.end() - n;
    auto u     = U.end() - 1;
    auto udiag = u;

    *bp-- /= *u;


    for (int i = n - 2; i >= 0; i--)
        {
        u     -= n;
        urow  -= n;
        udiag -= np1;

        int j = i + 1;

        double dotProduct = 0.0;
        for (auto uj = urow + j, end = urow + n, bj = b + j; uj < end; uj++, bj++)
            dotProduct += *uj * *bj;

        *bp = (*bp - dotProduct) / *udiag;
        --bp;
        }
}

void MathCache::forwardSubstitutionRow(DoubleMatrix& L, double* const b) {

    auto lrows  = L.getNumRows();
    auto lcols  = L.getNumCols();
    auto lcols1 = lcols+1;
    auto lrow   = L.begin();
    auto ldiag  = lrow;
    auto bp     = b;

    *bp++ /= *lrow;
    for (size_t i = 1; i < lrows; i++)
        {
        lrow  += lcols;
        ldiag += lcols1;

        double dotProduct = 0.0;
        for (auto lp = lrow, lend = lrow + i, bj = b; lp < lend; ++lp, ++bj)
            dotProduct += *lp * *bj;

        *bp = (*bp - dotProduct) / *ldiag;
        ++bp;
        }
}

void MathCache::computeLandU(DoubleMatrix& A, DoubleMatrix& L, DoubleMatrix& U) {

    auto n      = A.getNumRows();
    auto abegin = A.begin();
    auto end    = A.end();
    auto n1     = n + 1;
    auto ajj    = abegin;
    auto akj0   = abegin;

    for (size_t j = 0; j < n; j++)
        {
        auto akj = akj0;
        for (size_t k = 0; k < j; k++, akj += n)
            {
            auto k1  = k + 1;
            auto aij = A.getPointer(k1, j);
            auto aik = A.getPointer(k1, k);
            for (size_t i = k1; i < j; i++)
                {
                *aij -= *aik * *akj;
                aij  += n;
                aik  += n;
                }
            }

        akj = akj0;
        for (size_t k = 0; k < j; k++, akj += n)
            {   
            auto aij = ajj;
            auto aik = A.getPointer(j, k);
            for (size_t i = j; i < n; i++)
                {
                *aij -= *aik * *akj;
                aij  += n;
                aik  += n;
                }
        }

        auto amj = A.getPointer(j + 1, j);
        while (amj < end) 
            {
            *amj /= *ajj;
            amj  += n;
            }

        ajj += n1;
        ++akj0;
        }

    auto u = U.begin();
    auto l = L.begin();
    auto a = abegin;
    for (size_t row = 0; row < n; row++)
        {
        for (size_t col = 0; col < n; col++)
            {
            if (row <= col)
                {
                *u = *a;
                *l = row == col ? 1.0 : 0.0;
                }
            else
                {
                *l = *a;
                *u = 0.0;
                }

            ++l;
            ++u;
            ++a;
            }
        }
}

void MathCache::gaussianElimination(DoubleMatrix& A, DoubleMatrix& B, DoubleMatrix& X) {

    assert(A.getNumRows() == A.getNumCols());
    
    auto n    = A.getNumRows();
    auto b    = pushArray(n)->begin();
    auto bend = b + n;
    auto l    = pushMatrix(n);
    auto s    = pushMatrix(n);

    computeLandU(A, *l, *s);

    for (auto xrow = X.begin(), brow = B.begin(), xend = xrow + n; xrow < xend; ++xrow, ++brow)
        {
        for (auto bp = b, br = brow; bp < bend; ++bp, br += n)
            *bp = *br;

        /* Answer of Ly = b (which is solving for y) is copied into b. */
        forwardSubstitutionRow(*l, b);

        /* Answer of Ux = y (solving for x and the y was copied into b above)
           is also copied into b. */
        backSubstitutionRow(*s, b);

        for (auto bp = b, xr = xrow; bp < bend; ++bp, xr += n)
            *xr = *bp;
        }

    popMatrix(2);
    popArray();
}

void MathCache::multiply(DoubleMatrix& A, DoubleMatrix& B) {

    auto m = pushMatrix(A.getNumCols(), B.getNumRows());
    A.multiply(B, *m);
    B.copy(*m);
    popMatrix(1);
}

// Raise m to an integer power
void MathCache::power(DoubleMatrix& m, int power) {

    assert(m.getNumRows() == m.getNumCols());
    auto n = m.getNumRows();

    // Handle the common low-power cases directly
    switch (power) {
      case 0:
          m.setIdentity();
          break;

      case 1:
          break;

      case 2: 
          {
          auto t = pushMatrix(n);
          m.multiply(m, *t);
          m.copy(*t);
          popMatrix(1);
          }
          break;

      case 3: 
          {
          auto t = pushMatrix(n);
          m.multiply(m, *t);
          multiply(*t, m);
          popMatrix(1);
          }
          break;

      case 4: 
          {
          auto t = pushMatrix(n);
          m.multiply(m, *t);
          t->multiply(*t, m);
          popMatrix(1);
          }
          break;

      default: 
          {
          // As long as there are 1-bits left in the power integer
          // Calculate even powers of 2 
          // Accumulate the result in "m"

          auto t    = pushMatrix(n);
          auto even = pushMatrix(n);
          even->copy(m);
          bool first = false;

          for (;;)
              {
              if ((power & 1) > 0) 
                  {
                  if (first) 
                      {
                      m.multiply(*even, *t);
                      m.copy(*t);
                      }
                  else 
                      {
                      first = true;
                      m.copy(*even);
                      }
                  }
              power >>= 1;
              if (power > 0) 
                  {
                  even->multiply(*even, *t);
                  even->copy(*t);
                  }
              else
                  break;
              }

          popMatrix(2);
          break;

        }
    }
}

void MathCache::transpose(DoubleMatrix& m, DoubleMatrix& mT) {

    size_t n = m.getNumRows();
	for (size_t i=0; i<n; i++)
		for (size_t j=0; j<n; j++)
			mT(j,i) = m(i,j);
}

