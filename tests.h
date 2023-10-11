#ifndef BIO_TESTS_H_
#define BIO_TESTS_H_

#include "aa.h"
#include "cdn.h"
#include "dna.h"

#include "polymer.h"

namespace bio {
namespace test {

class test_failed_error : public std::runtime_error {
    public: test_failed_error(const std::string &what) : std::runtime_error(what) {}
};

void run_all();

void nts_from_string();
void cdns_from_string();
void aas_from_string();
void aas_from_nts();
void rc_nts();

};
};


#endif