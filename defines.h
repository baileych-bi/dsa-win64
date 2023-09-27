#ifndef DEFINES_H_
#define DEFINES_H_

#ifndef VERSION_STRING
#define VERSION_STRING "0.3.0"
#endif

#include <vector>

namespace bio {
    template<typename T>
    using vecvec = std::vector<std::vector<T>>;
};

#endif
