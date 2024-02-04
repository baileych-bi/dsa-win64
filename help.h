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

#ifndef BIO_HELP_H_
#define BIO_HELP_H_

#include <iomanip>
#include <iostream>

#include "aa.h"
#include "params.h"

#ifdef DSA_TARGET_WIN64
#include "local-getopt.h"
#elif defined(DSA_TARGET_LINUX)
#include <getopt.h>
#endif

namespace help {

using bio::Aas;

/** Help text for a particular command line option */
struct OptHelp {
    char cname = 0;               ///< The one-letter command line flag/
    const char *sname=nullptr;    ///< The "long option" name.
    const char *helptext=nullptr; ///< Help text/
};

/** Print a command line option and its help text with proper indentation */
void
print_option(std::ostream &os, const char *option, const char *helptext);

/** Print multiple command line options and their help text with proper indentation */
template<size_t N>
void
print_opthelp(const OptHelp (&help)[N]) {
    size_t maxw = 0;
    for (size_t i=0; i<N; ++i) {
        const OptHelp &h = help[i];
        size_t w = h.sname ? std::strlen(h.sname) : 0;
        if (w > maxw) maxw = w; 
    }
    maxw += 2;

    for (size_t i=0; i<N; ++i) {
        const OptHelp &h = help[i];
        std::cout << "  " << (h.cname ? "-" : " ") << (h.cname ? h.cname : ' ');
        std::cout << "  ";
        std::cout << (h.sname ? "--" : "  ");
        std::cout << std::setw(maxw) << std::left << (h.sname ? h.sname : " ");
        if (h.helptext) std::cout << h.helptext;
        std::cout << std::endl;
    }
}

/** Print the help text */
void
print_help();

/** Print the ascii codon table */
void
print_codon_table();

/** Use <getopt> to parse the command line options into Params */
Params
parse_argv(int argc, char **argv);

}; //namespace help
#endif