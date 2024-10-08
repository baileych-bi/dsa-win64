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

#include <unordered_map>

#include "params.h"

namespace help {

std::optional<CodonOutput>
codon_output_from_string(const char *s) {
    static const std::unordered_map<std::string, CodonOutput> lookup = {
        {"none",       CodonOutput::None},
        {"ascii",      CodonOutput::Ascii},
        {"horizontal", CodonOutput::Horizontal},
        {"vertical",   CodonOutput::Vertical}
    };

    std::optional<CodonOutput> co;
    std::string lc; for (; *s; ++s) lc.push_back(std::tolower(*s));
    auto ii = lookup.find(lc);
    if (ii != lookup.end()) co = ii->second;
    return co;
}

};