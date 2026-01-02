#include "Container.hpp"


BufferAllocator::BufferAllocator(void) {

    buffer = nullptr;
    currentSize = 0;
    maxSize = 0;
}

BufferAllocator::~BufferAllocator(void) {

    deallocate();
}

void BufferAllocator::deallocate(void) {

    if (buffer) 
        {
        free(buffer);
        buffer = nullptr;
        currentSize = 0;
        maxSize = 0;
    }
}

void BufferAllocator::allocate(size_t size) {

    if (size > maxSize)
        {
        deallocate();
        if (size > 0) 
            buffer = malloc(size);
        maxSize = size;
        }
    currentSize = size;
}

void BufferAllocator::setZero(void) {

    memset(buffer, 0, currentSize);
}

void BufferAllocator::copy(const BufferAllocator &b) {

    allocate(b.currentSize);
    memcpy(buffer, b.buffer, currentSize);
}

bool BufferAllocator::operator==(const BufferAllocator& b) const {

    return currentSize == b.currentSize && memcmp(buffer, b.buffer, currentSize) == 0;
}

bool BufferAllocator::operator!=(const BufferAllocator& b) const {

    return currentSize != b.currentSize ||memcmp(buffer, b.buffer, currentSize) != 0;
}


