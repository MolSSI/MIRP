/*! \file
 *
 * \brief Functions related to testing the Boys function
 */

#include "mirp_bin/boys_test.hpp"
#include "mirp_bin/test_common.hpp"
#include <mirp/kernels/boys.h>
#include <mirp/math.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <memory>

/* Anonymous namespace for some helper functions */
namespace {


/* Runs a Boys function test using interval arithmetic
 *
 * This just wraps the calculation of Boys function, increasing the
 * working precision until the desired target precision is reached. It
 * then handles the comparison with the reference data
 *
 * The number of failing tests is returned
 *
 * \todo This function is not exception safe
 */
long boys_run_test(const mirp::boys_data & data, long extra_m, long target_prec)
{
    using namespace mirp;

    long nfailed = 0;

    const int max_m = boys_max_m(data) + extra_m;

    arb_t t_mp, vref_mp;
    arb_init(t_mp);
    arb_init(vref_mp);

    arb_ptr F_mp = _arb_vec_init(max_m+1);

    for(const auto & ent : data.values)
    {
        /* 16 extra bits (~4-5 decimal digits) for safety */
        mirp_boys_target_prec_str(F_mp, ent.m + extra_m, ent.t.c_str(), target_prec+16);

        /* Round the reference value to the target precision */
        arb_set_str(vref_mp, ent.value.c_str(), target_prec);

        /* Rounding the reference value to the target precision results in
         * an interval. Does that interval contain our (more precise) result?
         */
        /* 1.) Test if the calculated value is within the error of the reference
         * 2.) If it's not, test if the calculated value is [0 +/ value] and is
         *     that zero within the target precision
         * 3.) If (2) is true, and the reference integral is exactly zero, then
         *     they are considered equal (and this block is not entered)
         */
        if(!arb_contains(vref_mp, F_mp + ent.m) &&
           !(arb_is_zero(vref_mp) && mirp_test_zero_prec(F_mp + ent.m, target_prec)))
        {
            std::cout << "Entry failed test: m = " << ent.m << " t = " << ent.t << "\n";
            char * s1 = arb_get_str(F_mp + ent.m, 1000, 0);
            char * s2 = arb_get_str(vref_mp, 1000, 0);
            std::cout << "   Calculated: " << s1 << "\n";
            std::cout << "    Reference: " << s2 << "\n\n";
            free(s1);
            free(s2);
            nfailed++;
        }
    }

    arb_clear(t_mp);
    arb_clear(vref_mp);
    _arb_vec_clear(F_mp, max_m+1);

    return nfailed;
}


/* Runs a Boys function test using double precision
 *
 * This just wraps the calculation of Boys function and the
 * comparison of the computed data with the reference values.
 *
 * The number of failing tests is returned
 */
long boys_run_test_d(const mirp::boys_data & data, long extra_m)
{
    using namespace mirp;

    long nfailed = 0;

    const int max_m = boys_max_m(data) + extra_m;
    std::vector<double> F_dbl(max_m+1);

    for(const auto & ent : data.values)
    {
        double t_dbl = std::strtod(ent.t.c_str(), nullptr);
        double vref_dbl = std::strtod(ent.value.c_str(), nullptr);

        mirp_boys_d(F_dbl.data(), ent.m+extra_m, t_dbl);

        if(vref_dbl != F_dbl[ent.m])
        {
            std::cout << "Entry failed test: m = " << ent.m << " t = " << ent.t << "\n";
            double reldiff = std::fabs(vref_dbl - F_dbl[ent.m]);
            reldiff /= std::fmax(std::fabs(vref_dbl), std::fabs(F_dbl[ent.m]));

            auto old_cout_prec = std::cout.precision(17);
            std::cout << "   Calculated: " << F_dbl[ent.m] << "\n";
            std::cout << "    Reference: " << vref_dbl << "\n";
            std::cout << "Relative Diff: " << reldiff << "\n\n";
            std::cout.precision(old_cout_prec);
            nfailed++;
        }
    }

    return nfailed;
}


/* Runs a Boys function test using 'exact' double precision
 *
 * This is just a simple test of the wrapper. The comparison
 * is not expected to be exact, since the inputs are in greater
 * than double precision. This results in rounding of the inputs,
 * and therefore the result as calculated in mirp_boys_exact
 * will differ from that calulcated purely in interval arithmetic.
 *
 * This, therefore, just ensures that the wrappers are written correctly.
 */
long boys_run_test_exact(const mirp::boys_data & data, long extra_m)
{
    using namespace mirp;
    using mirp::detail::almost_equal;

    long nfailed = 0;

    const int max_m = boys_max_m(data) + extra_m;
    std::vector<double> F_dbl(max_m+1);

    for(const auto & ent : data.values)
    {
        double t_dbl = std::strtod(ent.t.c_str(), nullptr);
        double vref_dbl = std::strtod(ent.value.c_str(), nullptr);

        mirp_boys_exact(F_dbl.data(), ent.m+extra_m, t_dbl);

        if(!almost_equal(vref_dbl, F_dbl[ent.m], 1e-13))
        {
            std::cout << "Entry failed test: m = " << ent.m << " t = " << ent.t << "\n";
            double reldiff = std::fabs(vref_dbl - F_dbl[ent.m]);
            reldiff /= std::fmax(std::fabs(vref_dbl), std::fabs(F_dbl[ent.m]));

            auto old_cout_prec = std::cout.precision(17);
            std::cout << "   Calculated: " << F_dbl[ent.m] << "\n";
            std::cout << "    Reference: " << vref_dbl << "\n";
            std::cout << "Relative Diff: " << reldiff << "\n\n";
            std::cout.precision(old_cout_prec);
            nfailed++;
        }
    }

    return nfailed;
}

} // close anonymous namespace


namespace mirp {

int boys_max_m(const boys_data & data)
{
    int max_m = 0;
    for(auto & ent : data.values)
        max_m = std::fmax(max_m, ent.m);
    return max_m;
}


boys_data boys_read_file(const std::string & filepath, bool is_input)
{
    using std::ifstream;

    ifstream infile(filepath, ifstream::in);
    if(!infile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for reading");

    std::string line;
    boys_data data;
    bool have_ndigits = false;

    while(std::getline(infile, line).good())
    {
        if(line.length() == 0)
            continue;
        else if(line[0] == '#')
            data.header += line + "\n";
        else
        {
            std::stringstream ss(line);
            ss.exceptions(std::stringstream::failbit);

            if(!is_input && !have_ndigits)
            {
                ss >> data.ndigits;
                have_ndigits = true;
            }
            else
            {
                boys_data_entry ent;
                ss >> ent.m >> ent.t;

                if(!is_input)
                    ss  >> ent.value;

                data.values.push_back(ent);
            }
        }
    }

    return data;
}


void boys_write_file(const std::string & filepath, const boys_data & data)
{
    using std::ofstream;

    ofstream outfile;
    outfile.open(filepath, ofstream::out | ofstream::trunc);

    if(!outfile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for writing");

    outfile.exceptions(ofstream::badbit | ofstream::failbit);

    outfile << data.header;
    outfile << data.ndigits << "\n";
    for(const auto & ent : data.values)
        outfile << ent.m << " " << ent.t << " " << ent.value << "\n";
}


long boys_run_test_main(const std::string & filepath,
                        const std::string & floattype,
                        long extra_m,
                        long target_prec)
{
    boys_data data;
    try {
        data = boys_read_file(filepath, false);
    }
    catch(std::exception & ex)
    {
        std::cout << "Unable to read data from file \"" << filepath << "\": ";
        std::cout << ex.what() << "\n";
        return 2;
    }

    std::cout << "Read " << data.values.size() << " values from " << filepath << "\n";

    long nfailed = 0;

    if(floattype == "interval")
        nfailed = boys_run_test(data, extra_m, target_prec);
    else if(floattype == "exact")
        nfailed = boys_run_test_exact(data, extra_m);
    else if(floattype == "double")
        nfailed = boys_run_test_d(data, extra_m);
    else
    {
        std::string err;
        err = "Unknown floating-point type \"" + floattype + "\"";
        throw std::runtime_error(err);
    }


    double percent_passed = 100.0 - (100.0 * (double)nfailed / (double)data.values.size());
    std::cout << nfailed << " / " << data.values.size() << " failed ("
              << percent_passed << "% passed)\n";

    return nfailed;
}


void boys_create_test(const std::string & input_filepath,
                      const std::string & output_filepath,
                      long ndigits, const std::string & header)
{
    boys_data data = boys_read_file(input_filepath, true);
    data.ndigits = ndigits;
    data.header += header;

    const int max_m = boys_max_m(data);

    arb_t t_mp;
    arb_ptr F_mp = _arb_vec_init(max_m+1);
    arb_init(t_mp);

    /* Target precision/accuracy, in bits, with a safety factor
       of 5 extra decimal digits */
    const slong target_prec = (ndigits+5) / MIRP_LOG_10_2;

    for(auto & ent : data.values)
    {
        mirp_boys_target_prec_str(F_mp, ent.m, ent.t.c_str(), target_prec);

        char * s = arb_get_str(F_mp + ent.m, ndigits, ARB_STR_NO_RADIUS);
        ent.value = s;
        free(s);
    }

    boys_write_file(output_filepath, data);
    arb_clear(t_mp);
    _arb_vec_clear(F_mp, max_m+1);
}


} // closing namespace mirp

