#include <iostream>
#include <limits>
#include <sstream>
#include "BitVector.hpp"
#include "Msg.hpp"

int BitVector::bitsPerUint = sizeof(unsigned) * 8;
unsigned BitVector::maxUInt = std::numeric_limits<unsigned>::max() ;
unsigned BitVector::highBit = (unsigned(1) << (BitVector::bitsPerUint-1));



/** BitVector constructor */
BitVector::BitVector(void) {

    numElements = 0;
    numUints = 0;
    mask = 0;
    v.clear();
}

/** BitVector copy constructor
* @param b reference to the BitVector to be copied */
BitVector::BitVector(const BitVector& b) {

    clone(b);
}

/** BitVector constructor
* @param n the number of elements in the bitset */
BitVector::BitVector(size_t n) {

    // calculate size of vector and allocate memory for bitfield
    numElements = n;
    numUints = (numElements / BitVector::bitsPerUint);
    if (numElements % BitVector::bitsPerUint != 0)
        numUints++;
    v.resize(numUints);
    for (size_t i=0; i<numUints; i++)
        v[i] = 0;
    
    // initialize the mask
    if (numElements % BitVector::bitsPerUint != 0)
        {
        mask = (BitVector::highBit >> (numElements % BitVector::bitsPerUint));
        mask--;
        mask ^= BitVector::maxUInt;
        }
    else
        {
        mask = BitVector::maxUInt;
        }
}

BitVector::BitVector(size_t n, bool def) {

    // calculate size of vector and allocate memory for bitfield
    numElements = n;
    numUints = (numElements / BitVector::bitsPerUint);
    if (numElements % BitVector::bitsPerUint != 0)
        numUints++;
    v.resize(numUints);
    for (size_t i=0; i<numUints; i++)
        v[i] = 0;
    
    // initialize the mask
    if (numElements % BitVector::bitsPerUint != 0)
        {
        mask = (BitVector::highBit >> (numElements % BitVector::bitsPerUint));
        mask--;
        mask ^= BitVector::maxUInt;
        }
    else
        {
        mask = BitVector::maxUInt;
        }

    // set all of the bits
    for (size_t i=0; i<numElements; i++)
        {
        if (def == true)
            set(i);
        else
            unset(i);
        }
}

/** BitVector subscript operator
* @param idx the bit to return */
bool BitVector::operator[](size_t idx) {

    if (idx >= numElements)
        Msg::error("Subscript index out-of-range");
    return (( v[idx/BitVector::bitsPerUint] & (BitVector::highBit >> (idx % BitVector::bitsPerUint))) != 0);
}

/** BitVector subscript operator
* @param idx the bit to return */
bool BitVector::operator[](size_t idx) const {

    if (idx >= numElements)
        Msg::error("Subscript index out-of-range");
    return (( v[idx/BitVector::bitsPerUint] & (BitVector::highBit >> (idx % BitVector::bitsPerUint))) != 0);
}

/** BitVector assignment operator
* @param b the bitset to equate this one to */
BitVector& BitVector::operator=(const BitVector& b) {

    if (this != &b)
        clone(b);
    return *this;
}

/** BitVector equality operator
* @param b the bitset to compare to this one */
bool BitVector::operator==(const BitVector& b) const {

    if (numElements != b.numElements)
        return false;
    for (int i=0; i<numUints; i++)
        if (v[i] != b.v[i])
            return false;
    return true;
}

/** BitVector not operator
* @param b the bitset whose bits are to be flipped */
BitVector operator!(const BitVector& b) {

    BitVector newB(b);
    newB.flip();
    return newB;
}

/** BitVector inequality operator
* @param lhs the left hand side bitset
* @param rhs the right hand side bitset */
bool operator!=(const BitVector& lhs, const BitVector& rhs) {

    return !(lhs == rhs);
}

/** BitVector less than or equal to operator
* @param b the bitset to compare to this one */
bool operator<=(const BitVector& a, const BitVector& b) {

    if (a > b)
        return false;
    else
        return true;
}

/** BitVector less than operator
* @param lhs the left hand side bitset
* @param rhs the right hand side bitset */
bool operator<(const BitVector& lhs, const BitVector& rhs) {

    if (lhs > rhs || lhs == rhs)
        return false;
    else
        return true;
}

/** BitVector greater than or equal to operator
* @param lhs the left hand side bitset
* @param rhs the right hand side bitset */
bool operator>=(const BitVector& lhs, const BitVector& rhs) {

    if (lhs > rhs || lhs == rhs)
        return true;
    else
        return false;
    
}

/** BitVector greater than operator
* @param b the right hand side bitset */
bool BitVector::operator>(const BitVector& b) const {

    size_t x = (numUints < b.numUints) ? numUints : b.numUints;
    for (size_t i=0; i<x; i++)
        {
        if (v[i] > b.v[i])
            return true;
        else if (v[i] < b.v[i])
            return false;
        }
    return false;
}

/** BitVector and operator
* @param a the right hand side bitset */
BitVector BitVector::operator&(const BitVector& a) const {

    if (a.numElements != numElements )
        return BitVector();
    else
        {
        BitVector b(numElements);
        for (size_t i=0; i<numUints; i++)
            b.v[i] = v[i] & a.v[i];
        return b;
        }
}

/** BitVector or operator
* @param a the right hand side bitset */
BitVector BitVector::operator|(const BitVector& a) const {

    if (a.numElements != numElements )
        return BitVector();
    else
        {
        BitVector b(numElements);
        for (size_t i=0; i<numUints; i++)
            b.v[i] = v[i] | a.v[i];
        return b;
        }
}

/** BitVector xor operator
* @param a the right hand side bitset */
BitVector BitVector::operator^(const BitVector& a) const {

    if (a.numElements != numElements )
        return BitVector();
    else
        {
        BitVector b(numElements);
        for (size_t i=0; i<numUints; i++)
            b.v[i] = v[i] ^ a.v[i];
        return b;
        }
}

/** Unary not */
BitVector& BitVector::operator~() {

    flip();
    return *this;
}

/** BitVector stream operator
* @param s the stream
* @param b the right hand side bitset */
std::ostream& operator<<(std::ostream& s, const BitVector& b) {

    s << "[" << b.bitString() << "]";
    return s;
}

/** BitVector and equals  operator
* @param lhs the left hand side bitset
* @param rhs the right hand side bitset */
BitVector& operator&=(BitVector& lhs, const BitVector& rhs) {

    if (lhs.size() == rhs.size())
        {
        for (size_t i=0; i<lhs.getNumUnints(); i++)
            {
            unsigned lvi = lhs.getUintElement(i);
            unsigned rvi = rhs.getUintElement(i);
            lhs.setUintElement(i, lvi &= rvi);
            }
        }
    return lhs;
}

/** BitVector or equals  operator
* @param lhs the left hand side bitset
* @param rhs the right hand side bitset */
BitVector& operator|=(BitVector& lhs, const BitVector& rhs) {

    if (lhs.size() == rhs.size())
        {
        for (size_t i=0; i<lhs.getNumUnints(); i++)
            {
            unsigned lvi = lhs.getUintElement(i);
            unsigned rvi = rhs.getUintElement(i);
            lhs.setUintElement(i, lvi |= rvi);
            }
        }
    return lhs;
}

/** BitVector xor equals  operator
* @param lhs the left hand side bitset
* @param rhs the right hand side bitset */
BitVector& operator^=(BitVector& lhs, const BitVector& rhs) {

    if (lhs.size() == rhs.size())
        {
        for (size_t i=0; i<lhs.getNumUnints(); i++)
            {
            unsigned lvi = lhs.getUintElement(i);
            unsigned rvi = rhs.getUintElement(i);
            lhs.setUintElement(i, lvi ^= rvi);
            }
        }
    return lhs;
}

/** Return a string with the bitset */
std::string BitVector::bitString(void) const {

    std::string str = "";
    for (size_t i=0; i<numElements; i++)
        {
        bool tf = (*this)[i];
        if (tf == true)
            str += "1";
        else
            str += "0";
        }
    return str;
}

/** Clone function
* @param b the BitVector to be copied */
void BitVector::clone(const BitVector& b) {

    numElements = b.numElements;
    numUints    = b.numUints;
    mask        = b.mask;
    v           = b.v;
}

/** Erase all bits */
void BitVector::clear(void) {

    numElements = 0;
    numUints = 0;
    mask = 0;
    v.clear();
}

/** Is the bitset of size zero */
bool BitVector::empty(void) const {

    return v.empty();
}

/** Flip the bit at the index
* @param idx the index of the bit to flip */
void BitVector::flip(size_t idx) {

    v[idx/BitVector::bitsPerUint] ^= ((BitVector::highBit >> (idx % BitVector::bitsPerUint)));
}

/** Flip all bits */
void BitVector::flip(void) {

    for (size_t i=0; i<numUints; i++)
        v[i] = v[i] ^ BitVector::maxUInt;
    v[numUints-1] &= mask;
}

/** The position of the first on bit
* @return The index of the first on bit */
size_t BitVector::getFirstSetBit(void) const {
    
    size_t i, j;
    for (i=0; i<numUints; i++)
        {
        if (v[i] != 0)
            break;
        }
    if (i == numUints)
        return -1;
    for (j=0; j<BitVector::bitsPerUint; j++)
        {
        if ((v[i] & (BitVector::highBit >> j)) != 0)
            break;
        }
    return (i * BitVector::bitsPerUint + j);

}

/** The number of on bits
* @return The number of on bits */
size_t BitVector::getNumberSetBits(void) const {
    
    size_t count = 0;
    for (size_t i=0; i<numElements; i++)
        {
        if ( (v[i/BitVector::bitsPerUint] & (BitVector::highBit >> (i % BitVector::bitsPerUint))) != 0 )
            count++;
        }
    return count;
}

/** Insert a bit
* @param pos the position in the bitset to insert
* @param tf the sign of the bit */
void BitVector::insert(size_t pos, bool tf) {

    // check position of insertion
    if (pos < 0 || pos > numElements)
        {
        std::cout << "Cannot insert bit at position " << pos << " because there are only " << numElements << " bits" << std::endl;
        return;
        }
    // resize, if necessary
    size_t n = numElements + 1;
    size_t nu = (n / BitVector::bitsPerUint);
    if (n % BitVector::bitsPerUint != 0)
        nu++;
    if (nu != numUints)
        v.push_back(0);
    
    numElements = n;
    numUints = nu;
    
    // initialize the mask
    if (numElements % BitVector::bitsPerUint != 0)
        {
        mask = (BitVector::highBit >> (numElements % BitVector::bitsPerUint));
        mask--;
        mask ^= BitVector::maxUInt;
        }
    else
        {
        mask = BitVector::maxUInt;
        }

    // set the bits
    for (size_t i=numElements-1; i>pos; i--)
        {
        if ( (*this)[i-1] == true )
            set(i);
        else
            unset(i);
        }
    if (tf == true)
        set(pos);
    else
        unset(pos);
}

/** The number of on bits
* @param idx the index of the bit
* @return  whether the bit is on */
bool BitVector::isSet(size_t idx) const {

    return ((v[idx/BitVector::bitsPerUint] & (BitVector::highBit >> (idx % BitVector::bitsPerUint))) != 0);
}

/** Print in a long format the bitset */
void BitVector::print(void) {

    std::string str = "";
    str += "BitVector\n";
    str += "   Size     = " + std::to_string(numElements) + "\n";
    str += "   Capacity = " + std::to_string(numUints*BitVector::bitsPerUint) + "\n   ";
    str += bitString();
    std::cout << str << std::endl;
}

/** Add a bit to the end of the set
* @param tf the value of the bit to add */
void BitVector::push_back(bool tf) {

    size_t n = numElements + 1;
    size_t nu = (n / BitVector::bitsPerUint);
    if (n % BitVector::bitsPerUint != 0)
        nu++;
    if (nu != numUints)
        v.push_back(0);
    
    numElements = n;
    numUints = nu;
    
    // initialize the mask
    if (numElements % BitVector::bitsPerUint != 0)
        {
        mask = (BitVector::highBit >> (numElements % BitVector::bitsPerUint));
        mask--;
        mask ^= BitVector::maxUInt;
        }
    else
        {
        mask = BitVector::maxUInt;
        }
    if (tf == true)
        set(numElements-1);
    else
        unset(numElements-1);
}

/** Add a bit to the front of the set
* @param tf the value of the bit to add */
void BitVector::push_front(bool tf) {

    size_t n = numElements + 1;
    size_t nu = (n / BitVector::bitsPerUint);
    if (n % BitVector::bitsPerUint != 0)
        nu++;
    if (nu != numUints)
        v.push_back(0);
    
    numElements = n;
    numUints = nu;
    
    // initialize the mask
    if (numElements % BitVector::bitsPerUint != 0)
        {
        mask = (BitVector::highBit >> (numElements % BitVector::bitsPerUint));
        mask--;
        mask ^= BitVector::maxUInt;
        }
    else
        {
        mask = BitVector::maxUInt;
        }

    for (size_t i=numElements-1; i>=1; i--)
        {
        if ( (*this)[i-1] == true )
            set(i);
        else
            unset(i);
        }
    if (tf == true)
        set(0);
    else
        unset(0);
}

/** Set all bist to on  */
void BitVector::set(void) {

    for (size_t i=0; i<numElements; i++)
        v[i/BitVector::bitsPerUint] |= (BitVector::highBit >> (i % BitVector::bitsPerUint));
}

/** Set a bit to on
* @param idx the index of the bit to turn on */
void BitVector::set(size_t idx) {

    if (idx >= numElements)
        Msg::error("Index out-of-range");
    v[idx/BitVector::bitsPerUint] |= (BitVector::highBit >> (idx % BitVector::bitsPerUint));
}

/** Set all bits to 0  */
void BitVector::unset(void) {

    for (size_t i=0; i<numElements; i++)
        v[i/BitVector::bitsPerUint] &= ((BitVector::highBit >> (i % BitVector::bitsPerUint))^BitVector::maxUInt);
}

/** Set a bit to 0
* @param idx the index of the bit to change */
void BitVector::unset(size_t idx) {

    if (idx >= numElements)
        Msg::error("Index out-of-range");
    v[idx/BitVector::bitsPerUint] &= ((BitVector::highBit >> (idx % BitVector::bitsPerUint))^BitVector::maxUInt);
}

/** Check the functionality of the class */
void BitVector::test(void) {

    // first check that edge cases are OK
    BitVector x1(10);
    x1.set(3);
    std::cout << x1 << " should be: (0001000000)" << std::endl;
    ~x1;
    std::cout << x1 << " should be: (1110111111)" << std::endl;
    BitVector x2(32);
    x2.set(31);
    std::cout << x2 << " should be: (00000000000000000000000000000001)" << std::endl;
    BitVector x3(33);
    x3.set(32);
    std::cout << x3 << " should be: (000000000000000000000000000000001)" << std::endl;
    
    // check operators
    BitVector x4(10);
    x4.set(4);
    std::cout << (x1 == x4) << " should be 0" << std::endl;
    BitVector x5(33);
    x5.set(32);
    std::cout << (x3 == x5) << " should be 1" << std::endl;
    std::cout << (x1 > x4) << " should be 1" << std::endl;
    std::cout << (x1 | x4) << " should be (0001100000)" << std::endl;
    std::cout << (x1 & x4) << " should be (0000000000)" << std::endl;
    std::cout << (x1 ^ x4) << " should be (0001100000)" << std::endl;
    std::cout << x1[3] << x2[31] << x3[32] << x4[4] << x5[32] << " should be 11111" << std::endl;
    std::cout << x1[0] << x2[14] << x3[30] << x4[8] << x5[29] << " should be 00000" << std::endl;
    
    BitVector x6(1);
    for (int i=0; i<100; i++)
        {
        x6.push_back(true);
        std::cout << x6 << std::endl;
        }

    BitVector x7(1);
    for (int i=0; i<100; i++)
        {
        x7.push_front(true);
        std::cout << x7 << std::endl;
        }
        
    BitVector x8(32);
    x8.insert(1, true);
    std::cout << x8 << std::endl;
    x8.print();

    BitVector x9(32);
    x9.insert(33, true);
    std::cout << x9 << std::endl;
    x9.print();
}
