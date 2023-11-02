#ifndef Container_H
#define Container_H

#include <iomanip>
#include <iostream>
#include <string>
#include <assert.h>

#define ForElements(v) for (auto v = this->begin(), end = this->end(); v < end; ++v)
#define ForLeftRight(a) const auto* right = a.begin(); for (auto left = this->begin(), end = this->end(); left < end; ++left, ++right)

#ifndef _CONSOLE
#pragma mark - BufferTemplate Definition -
#endif

class BufferAllocator {

    public:
        explicit BufferAllocator(void);
        virtual ~BufferAllocator(void);
        void setZero(void);
        bool operator==(const BufferAllocator &b) const;
        bool operator!=(const BufferAllocator &b) const;

    protected:
        virtual void deallocate();
        void allocate(size_t size);
        void copy(const BufferAllocator& b);
        const void* getBuffer() const {return buffer;}
        void* getBuffer() { return buffer; }

    private:
        void*  buffer;
        size_t maxSize,
               currentSize;
};


template<typename T>
class BufferTemplate: public BufferAllocator {

    public:
        BufferTemplate(void) {
        
            endBuffer = NULL;
            numElements = 0;
        }

        explicit BufferTemplate(const BufferTemplate<T>& a) {
        
            copy(a);
        }

        explicit BufferTemplate(size_t elements) {
        
            create(elements);
        }

        void create(size_t elements) {
        
            allocate(elements * sizeof(T));
            numElements = elements;
            endBuffer   = begin() + elements;
        }

        void copy(const BufferTemplate<T>& b) {
        
            BufferAllocator::copy(b);
            numElements = b.numElements;
            endBuffer   = begin() + numElements;
        }

        void fill(T value) {
            ForElements(e)
                *e = value;
        }

        void add(T c) {

            ForElements(e)
                *e += c;
        }

        void add(const BufferTemplate<T>& a) {

            ForLeftRight(a)
                *left += *right;
        }

        void multiply(T scalar) {

            ForElements(e)
                *e *= scalar;
        }

        void multiply(T scalar, BufferTemplate<T>& result) const {
        
            result.create(numElements);
            auto r = result.begin();
            ForElements(e)
                *r++ = *e * scalar;
        }

        void operator-=(T scalar) {

            ForElements(e)
                *e -= scalar;
        }

        void operator/=(T scalar) {

            ForElements(e)
                *e /= scalar;
        }

        void operator|=(T scalar) {

            ForElements(e)
                *e |= scalar;
        }

        void operator^=(T scalar) {

            ForElements(e)
                *e ^= scalar;
        }

        void operator-=(const BufferTemplate<T>& a) {

            ForLeftRight(a)
                *left -= *right;
        }

        const T*    begin(void) const { return (T*)this->getBuffer(); }
        T*          begin(void) { return (T*)this->getBuffer(); }
        const T*    end(void) const { return (T*)endBuffer; }
        T*          end(void) { return (T*)endBuffer; }
        size_t      size(void) { return numElements; }

        T*          endBuffer;
        size_t      numElements;
};



#ifndef _CONSOLE
#pragma mark - ArrayTemplate Definition -
#endif

template<typename T>
class ArrayTemplate: public BufferTemplate<T> {

    public:
        ArrayTemplate(void) {
        
        }

        ArrayTemplate(size_t elements) : BufferTemplate<T>(elements) {
        
        }

        ArrayTemplate(const ArrayTemplate<T>& a) : BufferTemplate<T>(a) {
        
        }

        ArrayTemplate<T>& operator=(const ArrayTemplate<T>& rhs) {
        
            if (this != &rhs)
                BufferTemplate<T>::copy(rhs);
            return *this;
        }

        T operator[](size_t i) const {
        
            return this->begin()[i]; 
        }

        T getValue(size_t i) const { return this->begin()[i]; }

        void copy(const ArrayTemplate<T>& other) {
        
            BufferTemplate<T>::copy(other);
        }
};



#ifndef _CONSOLE
#pragma mark - MatrixTemplate Definition -
#endif

template<typename T>
class MatrixTemplate : public BufferTemplate<T> {

    public:
        explicit MatrixTemplate(void) {
        
            numRows = 0;
            numCols = 0;
        }

        explicit MatrixTemplate(size_t nr, size_t nc) : BufferTemplate<T>::BufferTemplate(nr * nc) {
        
            numRows = nr;
            numCols = nc;
            this->setZero();
        }

        explicit MatrixTemplate(const MatrixTemplate<T>& m) : BufferTemplate<T>::BufferTemplate(m) {

            numRows = m.numRows;
            numCols = m.numCols;
        }

        MatrixTemplate<T>& operator=(const MatrixTemplate<T>& rhs) {
        
            if (this != &rhs)
                copy(rhs);
            return *this;
        }

        T& operator()(size_t r, size_t c) {
        
            return this->begin()[r * numCols + c];
        }

        const T& operator()(size_t r, size_t c) const {
        
            return this->begin()[r * numCols + c];
        }

        size_t getNumRows() const {return numRows;}
        size_t getNumCols() const { return numCols; }
        T* getRow(int r) { return this->begin() + r * numCols; }

        T* getPointer(size_t r, size_t c) {

            return this->begin() + (r * numCols + c);
        }

        T getValue(size_t r, size_t c) const {

            return this->begin()[r * numCols + c];
        }

        void setValue(size_t r, size_t c, T value) {

            this->begin()[r * numCols + c] = value;
        }

        void print(void) {

            std::cout << "{";
            for (int i = 0; i < numRows; i++)
                {
                std::cout << "{";
                for (int j = 0; j < numCols; j++)
                    {
                    std::cout << getValue(i, j);
                    if (j != numCols - 1)
                        std::cout << ",";
                    }
                if (i != numRows - 1)
                    std::cout << "}," << std::endl;
                }
            std::cout << "}}" << std::endl;
            }

        void print(std::string title) {
        
            std::cout << title << std::endl;
            print();
        }

        void copy(const MatrixTemplate<T>& other) {
        
            BufferTemplate<T>::copy(other);
            numRows = other.numRows;
            numCols = other.numCols;
        }

        void create(size_t nr, size_t nc) {
        
            if (numRows != nr || numCols != nc) 
                { 
                BufferTemplate<T>::create(nr * nc);
                numRows = nr;
                numCols = nc;
                }
        }

        void create(size_t n) {
        
            if (numRows != n || numCols != n) 
                {
                BufferTemplate<T>::create(n * n);
                numRows = n;
                numCols = n;
                }
        }

        void setIdentity(void) {

            assert(numRows == numCols);
            this->setZero();
            auto cols1 = getNumCols() + 1;
            for (auto end = this->end(), c = this->begin(); c < end; c += cols1)
                *c = 1;
        }

        double maxDiagonal(void) {
        
            assert(numRows == numCols);
            auto cols1 = getNumCols() + 1;
            T max = 0;
            for (auto end = this->end(), c = this->begin(); c < end; c += cols1)
                {
                if (*c > max)
                    max = *c;
                }
            return max;
        }

        bool operator==(const MatrixTemplate<T>& m) const {

            if (numRows != m.numRows || numCols != m.numCols)
                return false;
            return BufferTemplate<T>::operator==(m);
        }

        bool operator!=(const MatrixTemplate<T>& m) const {

            if (numRows != m.numRows || numCols != m.numCols)
                return true;
            return BufferTemplate<T>::operator!=(m);
        }

        void transpose(MatrixTemplate<T>& result) {
        
            assert(result.numRows == numCols && result.numCols == numRows);
            auto cell = this->begin();
            for (size_t r = 0; r < numRows; ++r)
                {
                for (int c = 0; c < numCols; ++c)
                    result.setValue(c, r, *cell++);
                }
        }

        void add(const MatrixTemplate<T>& m) {

            assert(numRows == m.numRows && numCols == m.numCols);
            BufferTemplate<T>::add(m);
        }

        void add(const MatrixTemplate<T>& m, MatrixTemplate<T>& result) const {

            assert(numRows == m.numRows && numCols == m.numCols);
            BufferTemplate<T>::add(m, result);
        }

        void multiply(T scalar) {
        
            BufferTemplate<T>::multiply(scalar);
        }

        void multiply(T scalar, MatrixTemplate<T>& result) const {

            assert(result.numRows == numRows && result.numCols == numCols);
            auto r = result.begin();
            ForElements(e)
                * r++ = *e * scalar;
        }

        void multiply(const MatrixTemplate<T>& m, MatrixTemplate<T>& result) const {

            assert(numCols == m.numRows);
            assert(result.numRows == numCols && result.numCols == numRows);

            auto cols = getNumCols();
            auto mcols = m.getNumCols();
            auto row = this->begin();
            auto end = this->end();
            auto r = result.begin();
            auto mbegin = m.begin();

            while (row < end)
                {
                auto mrow = mbegin;
                auto mend = mbegin + mcols;
                while (mrow < mend)
                    {
                    T sum = 0;
                    for (auto col = row, cend = row + cols, mcol = mrow; col < cend; ++col, mcol += mcols) 
                        sum += *col * *mcol;

                    *r++ = sum;
                    ++mrow;
                    }

                row += cols;
                }
        }

        void divideByPowerOfTwo(int power) {
        
            if (power > 0)
                {
                auto factor = 1L << power;
                multiply(1.0 / (double)factor);
                }
        }

    protected:
        size_t numRows,
               numCols;
};



#ifndef _CONSOLE
#pragma mark - Type Definition -
#endif


typedef ArrayTemplate<int>    IntArray;
typedef ArrayTemplate<double> DoubleArray;

typedef MatrixTemplate<int>    IntMatrix;
typedef MatrixTemplate<double> RealMatrix;

#endif
