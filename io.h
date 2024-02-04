/*
Copyright 2024, The Broad Institute of MIT and Harvard

Original Author: Charles C Bailey

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef CCB_IO_H_
#define CCB_IO_H_

#include <cassert>
#include <exception>
#include <memory>
#include <filesystem>

namespace bio {

namespace fs = std::filesystem;

/** The stage of the file mapping process where an exception occurred. */
enum class MappingStage {
    Open,
    Stat,
    Map
};

/** Error information for a failed file map event. */
class MappingException : public std::exception {
public:
    MappingException(MappingStage stage, int error)
        : std::exception()
        , stage(stage)
        , error(error) {}

    MappingStage stage = MappingStage::Open; ///< The stage of the mapping process where the error occurred
    int error = 0;                           ///< The errno returned by open(), stat(), or mmap()
};

struct ConstMappingImplementation; //hide platform-specific mmap behavior in here

/** Efficiently get the contents of a file as const char [].
  * Uses open(), stat(), and mmap() for Linux or the MapViewOfFile et al. APIs for Windows
  */
struct ConstMapping {
    ConstMapping();
    ConstMapping(const ConstMapping &)=delete;            ///< Mappings should be unique per file.
    ConstMapping &operator=(const ConstMapping &)=delete; ///< Mappings should be unique per file.
    ConstMapping(ConstMapping &&m) noexcept;
    ConstMapping &operator=(ConstMapping &&m) noexcept;

    ~ConstMapping();

    static ConstMapping map(const fs::path& path); ///< Map the given file to memory.

    const char* c_str() const;
    char operator[](size_t i) const;
    size_t size() const;

    const char* begin() const;
    const char* end()   const;
    void unmap();

private:
    ConstMapping(std::unique_ptr<ConstMappingImplementation>&& impl);
    std::unique_ptr<ConstMappingImplementation> impl_;
};

/** Advance pointer into buffered file over given number of newlines.
  *
  * @param cur the starting pointer in the buffer.
  * @param n the number of newlines to skip.
  * @param end pointer to last+1 character of the file.
  *
  * @return a pointer to the character after the nth newline, or end if end of buffer reached first.
 */
const char *
next_lines(const char *cur, size_t n, const char *end);

/** Advance arbitrary pointer into buffered fastq file to start of next record.
  *
  * @param cur the starting pointer in the buffer.
  * @param begin pointer to first character in buffer.
  * @param end pointer to last+1 character of the file.
  *
  * @return a pointer to the first character of the next record, or end if end of buffer reached first.
 */
const char *
seek_next(const char *cur, const char *begin, const char *end);

}; //namespace bio

#endif
