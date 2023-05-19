#include "io.h"

#include <iostream>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef DSA_TARGET_WIN64
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef min //undefine infuriating min and max macros from windows.h
#undef max //undefine infuriating min and max macros from windows.h

namespace bio {
struct ConstMappingImplementation {
    HANDLE fd_ = INVALID_HANDLE_VALUE;
    HANDLE mapping_ = INVALID_HANDLE_VALUE;
    const char* buf_ = nullptr;
    size_t size_ = 0;

    ConstMappingImplementation(const std::filesystem::path &path) {
        try {
            fd_ = CreateFile(
                path.c_str(),
                GENERIC_READ,
                0,
                NULL,
                OPEN_EXISTING,
                FILE_ATTRIBUTE_NORMAL,
                NULL
            );

            if (fd_ == INVALID_HANDLE_VALUE) {
                throw MappingException(MappingStage::Open, errno);
            }

            LARGE_INTEGER file_size;

            if (!GetFileSizeEx((HANDLE)fd_, &file_size)) {
                throw MappingException(MappingStage::Stat, errno);
            }

            size_ = file_size.QuadPart;

            mapping_ = CreateFileMapping(
                (HANDLE)fd_,
                NULL,
                PAGE_READONLY,
                0,
                0,
                NULL
            );

            if (mapping_ == INVALID_HANDLE_VALUE) {
                throw MappingException(MappingStage::Map, errno);
            }

            buf_ = reinterpret_cast<const char *>(MapViewOfFile(
                mapping_,
                FILE_MAP_READ,
                0,
                0,
                0
            ));

            if (!buf_) {
                throw MappingException(MappingStage::Map, errno);
            }
        }
        catch (const MappingException& ex) {
            unmap();
            throw ex;
        }
    }

    void unmap() {
        if (buf_) UnmapViewOfFile(buf_);

        buf_ = nullptr;
        size_ = 0;

        if (mapping_ != INVALID_HANDLE_VALUE) {
            CloseHandle(mapping_);
            mapping_ = INVALID_HANDLE_VALUE;
        }
        if (fd_ != INVALID_HANDLE_VALUE) {
            CloseHandle(fd_);
            fd_ = INVALID_HANDLE_VALUE;
        }
    }

    const char* c_str() const { return buf_; }
    char operator[](size_t i) const { return *(buf_ + i); }
    size_t size() const { return size_; }

    ~ConstMappingImplementation() { unmap(); }
};
}; //namespace bio

#elif defined(DSA_TARGET_LINUX) //#ifdef DSA_TARGET_WIN64_
#include <sys/mman.h>
#include <unistd.h>

namespace bio {
struct ConstMappingImplementation {
    const char* buf_ = nullptr;
    size_t size_ = 0;
    int fd_ = -1;

    ConstMappingImplementation(const std::filesystem::path& path) {
        try {
            fd_ = open(path.c_str(), O_RDONLY);

            if (fd_ < 0) throw MappingException(MappingStage::Open, errno);

            struct stat statbuf;
            int err = fstat(fd_, &statbuf);
            if (err < 0) throw MappingException(MappingStage::Stat, errno);

            size_ = statbuf.st_size;
            buf_ = static_cast<const char*>(mmap(NULL, statbuf.st_size, PROT_READ, MAP_PRIVATE, fd_, 0));
            if (MAP_FAILED == buf_) throw MappingException(MappingStage::Map, errno);

            close(fd_);
        }
        catch (const MappingException& ex) {
            unmap();
            throw ex;
        }
    }

    ~ConstMappingImplementation() { unmap(); }

    const char* c_str() const { return buf_; }

    char operator[](size_t i) const { return *(buf_ + i); }

    size_t size() const { return size_; }

    void unmap() {
        if (buf_) munmap(reinterpret_cast<void*>(const_cast<char*>(buf_)), size_);
        buf_ = nullptr;
        size_ = 0;

        if (fd_ >= 0) close(fd_);
        fd_ = -1;
    };
};
}; //namespace bio

#endif //#elif defined(DSA_TARGET_LINUX_)

namespace bio {

ConstMapping
ConstMapping::map(const fs::path& path) {
    return ConstMapping(std::make_unique<ConstMappingImplementation>(path));
} ///< Map the given file to memory.

ConstMapping::ConstMapping() : impl_(nullptr) {}

ConstMapping::ConstMapping(std::unique_ptr<ConstMappingImplementation>&& impl)
    : impl_(std::move(impl)) {}

ConstMapping::~ConstMapping() { unmap(); }

ConstMapping::ConstMapping(ConstMapping&& rhs) noexcept
    : impl_(std::move(rhs.impl_)) {}

ConstMapping &
ConstMapping::operator=(ConstMapping&& rhs) noexcept {
    std::swap(impl_, rhs.impl_);
    return *this;
}

const char* ConstMapping::c_str() const { return impl_->c_str(); }
char ConstMapping::operator[](size_t i) const { return (*impl_)[i]; }
size_t ConstMapping::size() const { return impl_->size(); }

const char* ConstMapping::begin() const { return impl_->c_str(); }
const char* ConstMapping::end()   const { return begin() ? begin() + size() : begin(); }

void ConstMapping::unmap() { impl_->unmap(); }

const char *
next_lines(const char *cur, size_t n, const char *end) {
    size_t newline_count = 0;
    for (; cur != end; ++cur) {
        if (*cur == '\n') {
            if (newline_count == n) return cur + 1;
            ++newline_count;
        }
    }
    return end;
}

const char *
seek_next(const char *cur, const char *begin, const char *end) {
    for (; cur != end; ++cur) {
        if (*cur == '+') {
            if (cur+1 == end) return end;
            if (cur != begin && *(cur-1) == '\n' && *(cur+1) == '\n') return next_lines(cur, 1, end);
        }
    }
    return end;
}

};