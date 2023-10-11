#include "tests.h"

namespace bio {
namespace test {

void
run_all() {
    nts_from_string();
    cdns_from_string();
    aas_from_string();
    aas_from_nts();
    rc_nts();
}

void
nts_from_string() {
     const char8_t *utf8 = u8"Lorem ipsum dolor sit amet, consectetur adipiscing elit. Donec tincidunt, augue nec mattis porta,"
        u8"The quick brown fox jumped over the lazy dog"
        u8" !\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
        u8"⑅⬄―✦⸷➞a⥼⡹∎♒t┩∡⾯⚯⌭c⛶ⅳⰯg⃶♼⿜➞⬺⥱✽⏠ⲵ⎎⧒⍴⠰ⴛ⣋ⱗ⥝▿ⷅⷶ⪭✢⚇⮢☔⻊╣↟✏ⓟ␹⪍⺼⫏⤮₷△⌬⭓➊⓰⾧⿟␅⹳⎄Ⱁ⯾⯴✠⹸⎵♨ⴐA⚍⣈≭Cⰽ≱⬉TⰗⱻG⍸⎑⡙↉∨⠉₟⿁Ⲕ⏊ⷮ⇔⦝⋻⚹⸭⁆⤤⭋⑜*";
    const std::string input = reinterpret_cast<const char *>(utf8);
    const std::string expected = "TATCNCTTACNGTNCTNCNTAGNCATTTATCNTAGACGNTACGNTATCGACTG";

    {
        Nts nts(input);
        if (nts.as_string_view() != expected) throw test_failed_error("Nts(const std::string &) failed");
    }

    {
        std::string s = input;
        Nts nts(std::move(s));
        if (nts.as_string_view() != expected) throw test_failed_error("Nts(std::string &&) failed");
        
        nts.clear();
        if (nts.as_string_view() != "") throw test_failed_error("Nts::clear() failed");

        nts = input;
        if (nts.as_string_view() != expected) throw test_failed_error("Nts::operator=(const std::string &) failed");

        s = input;
        nts = std::move(s);
        if (nts.as_string_view() != expected) throw test_failed_error("Nts::operator+=(std::string &&) failed");
    }
}

void
aas_from_string() {
    const char8_t *utf8 = u8"Lorem ipsum dolor sit amet, consectetur adipiscing elit. Donec tincidunt, augue nec mattis porta,"
        u8"The quick brown fox jumped over the lazy dog"
        u8" !\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
        u8"⑅⬄―✦⸷➞a⥼⡹∎♒t┩∡⾯⚯⌭c⛶ⅳⰯg⃶♼⿜➞⬺⥱✽⏠ⲵ⎎⧒⍴⠰ⴛ⣋ⱗ⥝▿ⷅⷶ⪭✢⚇⮢☔⻊╣↟✏ⓟ␹⪍⺼⫏⤮₷△⌬⭓➊⓰⾧⿟␅⹳⎄Ⱁ⯾⯴✠⹸⎵♨ⴐA⚍⣈≭Cⰽ≱⬉TⰗⱻG⍸⎑⡙↉∨⠉₟⿁Ⲕ⏊ⷮ⇔⦝⋻⚹⸭⁆⤤⭋⑜*";
    
    const std::string input = reinterpret_cast<const char *>(utf8);
    const std::string expected = "LREMIPSMDLRSITAMETCNSECTETRADIPISCINGELITDNECTINCIDNTAGENECMATTISPRTATHEQICKRWNFMPEDVERTHELAYDG*ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYATCGACTG*";

    {
        Aas aas(input);
        if ( aas.as_string_view() != expected ) throw test_failed_error("Aas(const std::string &) failed : result was '" + std::string(aas.as_string_view()) + "'");
    }

    {
        std::string s = input;
        Aas aas(std::move(s));
        if ( aas.as_string_view() != expected ) throw test_failed_error("Aas(std::string &&) failed");

        aas.clear();
        if ( aas.as_string_view() != "" ) throw test_failed_error("Aas::clear() failed");

        aas = input;
        if ( aas.as_string_view() != expected ) throw test_failed_error("Aas::operator=(const std::string &) failed");

        s = input;
        aas = std::move(s);
        if ( aas.as_string_view() != expected ) throw test_failed_error("Aas::operator+=(std::string &&) failed");
    }
}

void aas_from_nts() {
    Nts nts = "AAAAACAATAAGACAACCACTACGATAATCATTATGAGAAGCAGTAGGCAACACCATCAGCCACCCCCTCCGCTACTCCTTCTGCGACGCCGTCGGTAATACTATTAGTCATCCTCTTCGTTATTCTTTTTGTGATGCTGTTGGGAAGACGATGAGGCAGCCGCTGCGGTAGTCGTTGTGGGAGGCGGTGGG";
    const std::string expected = "KNNKTTTTIIIMRSSRQHHQPPPPLLLLRRRR*YY*SSSSLFFL*CCWEDDEAAAAVVVVGGGG";

    {
        Nts temp = nts;
        Aas aas = std::move(temp);

        if (aas.as_string_view() != expected) throw test_failed_error("Aas(Nts &&) failed");
    }

    {
        Nts temp = nts;
        temp.exo(3, 3);

        Aas aas = std::move(temp);

        if (aas.as_string_view() != expected.substr(1, expected.size() - 2)) throw test_failed_error("Aas(Nts &&) failed after Nts::exo() called");
    }
}

void rc_nts() {
    const Nts fw = "TNCAANNCTCNNCGAGGNCAGNTCNACTAGGTGCTNACCGGTGNCAAAACTNTCNTGTNNGCCNAGAAGNCCTATNGCGAANGTGATCGCTGNNTTTAAT";
    std::string expected = "ATTAAANNCAGCGATCACNTTCGCNATAGGNCTTCTNGGCNNACANGANAGTTTTGNCACCGGTNAGCACCTAGTNGANCTGNCCTCGNNGAGNNTTGNA";

    {
        Nts temp = fw;
        temp.reverse_complement();
        if (temp.as_string_view() != expected) throw test_failed_error("Nts::reverse_complement() failed");
    }

    {
        Nts temp = fw;
        temp.exo(1, 1);
        temp.reverse_complement();
        if (temp.as_string_view() != expected.substr(1, expected.size() - 2)) throw test_failed_error("Nts::reverse_complement() failed after Nts::exo() called");
    }
}

void
cdns_from_string() {
    const char8_t *utf8 = u8"Lorem ipsum dolor sit amet, consectetur adipiscing elit. Donec tincidunt, augue nec mattis porta,"
        u8"The quick brown fox jumped over the lazy dog"
        u8" !\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
        u8"⑅⬄―✦⸷➞a⥼⡹∎♒t┩∡⾯⚯⌭c⛶ⅳⰯg⃶♼⿜➞⬺⥱✽⏠ⲵ⎎⧒⍴⠰ⴛ⣋ⱗ⥝▿ⷅⷶ⪭✢⚇⮢☔⻊╣↟✏ⓟ␹⪍⺼⫏⤮₷△⌬⭓➊⓰⾧⿟␅⹳⎄Ⱁ⯾⯴✠⹸⎵♨ⴐA⚍⣈≭Cⰽ≱⬉TⰗⱻG⍸⎑⡙↉∨⠉₟⿁Ⲕ⏊ⷮ⇔⦝⋻⚹⸭⁆⤤⭋⑜*";

    const std::string input = reinterpret_cast<const char *>(utf8);
    std::string expected;
    for (char c : input) if (Cdn::BIAS <= c && c < Cdn::BIAS + 64) expected.push_back(c);

    {
        Cdns cdns(input);
        if (cdns.as_string_view() != expected) throw test_failed_error("Cdns(const std::string &) failed");
    }

    {
        std::string s = input;
        Cdns cdns(std::move(s));
        if (cdns.as_string_view() != expected) throw test_failed_error("Cdns(std::string &&) failed");

        cdns.clear();
        if (cdns.as_string_view() != "") throw test_failed_error("Cdns::clear() failed");

        cdns = input;
        if (cdns.as_string_view() != expected) throw test_failed_error("Cdns::operator=(const std::string &) failed");

        s = input;
        cdns = std::move(s);
        if (cdns.as_string_view() != expected) throw test_failed_error("Cdns::operator+=(std::string &&) failed");
    }
}

}; //namespace test
}; //namespace bio
