#ifndef BitVector_hpp
#define BitVector_hpp

#include <string>
#include <vector>


class BitVector {

    public:
                                BitVector(void);
                                BitVector(const BitVector& b);
                                BitVector(size_t n);
                                BitVector(size_t n, bool def);
        bool                    operator[](size_t idx);
        bool                    operator[](size_t idx) const;
        BitVector&              operator=(const BitVector& b);
        bool                    operator==(const BitVector& b) const;
        bool                    operator>(const BitVector& b) const;
        BitVector               operator&(const BitVector& b) const;
        BitVector               operator|(const BitVector& b) const;
        BitVector               operator^(const BitVector& b) const;
        BitVector&              operator~();
        friend std::ostream&    operator<<(std::ostream& s, const BitVector& b);
        std::string             bitString(void) const;
        void                    clear(void);
        size_t                  getFirstSetBit(void) const;
        bool                    empty(void) const;
        void                    flip(size_t idx);
        void                    flip(void);
        size_t                  getNumberSetBits(void) const;
        size_t                  getNumUnints(void) const { return numUints; }
        unsigned                getUintElement(size_t idx) const { return v[idx]; }
        void                    insert(size_t pos, bool tf);
        bool                    isSet(size_t idx) const;
        void                    print(void);
        void                    push_back(bool tf);
        void                    push_front(bool tf);
        void                    set(void);
        void                    set(size_t idx);
        void                    setUintElement(size_t idx, unsigned x) { v[idx] = x; }
        size_t                  size(void) const { return numElements; }
        void                    test(void);
        void                    unset(void);
        void                    unset(size_t idx);

    private:
        void                    clone(const BitVector& b);
        static int              bitsPerUint;
        static unsigned         maxUInt;
        static unsigned         highBit;
        size_t                  numElements;
        size_t                  numUints;
        unsigned                mask;
        std::vector<unsigned>   v;
};

        BitVector               operator!(const BitVector& b);
        BitVector&              operator&=(BitVector& lhs, const BitVector& rhs);
        BitVector&              operator|=(BitVector& lhs, const BitVector& rhs);
        BitVector&              operator^=(BitVector& lhs, const BitVector& rhs);
        bool                    operator!=(const BitVector& lhs, const BitVector& rhs);
        bool                    operator< (const BitVector& lhs, const BitVector& rhs);
        bool                    operator<=(const BitVector& lhs, const BitVector& rhs);
        bool                    operator>=(const BitVector& lhs, const BitVector& rhs);

#endif
