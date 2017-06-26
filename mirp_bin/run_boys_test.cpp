
#include "test_files.hpp"
#include <mirp/mirp.h>
#include <memory>
#include <cstdlib>
#include <iostream>

using namespace mirp;

static long run_boys_test_mp(const boys_data & data, long target_prec, long working_prec)
{
    int maxm = 0;
    for(auto & it : data.values)
        maxm = std::max(maxm, it.m);

    mpfr_t t_mp, vref_mp;
    mpfr_init2(t_mp, working_prec);
    mpfr_init2(vref_mp, target_prec);

    auto F_mp = std::unique_ptr<mpfr_t[]>(new mpfr_t[maxm+1]);
    for(int i = 0; i <= maxm; i++)
        mpfr_init2(F_mp[i], target_prec);

    long nfailed = 0;

    for(const auto & it : data.values)
    {
        mpfr_set_str(t_mp, it.t.c_str(), 10, MPFR_RNDN);
        mpfr_set_str(vref_mp, it.value.c_str(), 10, MPFR_RNDN);
        mirp_boys_mp(F_mp.get(), it.m, t_mp, working_prec);

        if(!mpfr_equal_p(F_mp[it.m], vref_mp))
        {
            std::cout << "Entry failed test: m = " << it.m << " t = " << it.t << "\n";
            mpfr_printf("  Calculated: %Re\n", F_mp[it.m]);
            mpfr_printf("   Reference: %Re\n", vref_mp);
            nfailed++;
        }
    }

    mpfr_clear(t_mp);
    mpfr_clear(vref_mp);
    for(int i = 0; i <= maxm; i++)
        mpfr_clear(F_mp[i]);

    return nfailed;
}

static long run_boys_test_interval(const boys_data & data, long target_prec, long working_prec)
{
    int maxm = 0;
    for(auto & it : data.values)
        maxm = std::max(maxm, it.m);

    arb_t t_mp, vref_mp;
    arb_init(t_mp);
    arb_init(vref_mp);

    auto F_mp = std::unique_ptr<arb_t[]>(new arb_t[maxm+1]);
    for(int i = 0; i <= maxm; i++)
        arb_init(F_mp[i]);

    long nfailed = 0;

    for(const auto & it : data.values)
    {
        arb_set_str(t_mp, it.t.c_str(), working_prec);
        arb_set_str(vref_mp, it.value.c_str(), target_prec);
        mirp_boys_interval(F_mp.get(), it.m, t_mp, working_prec);

        arb_set_round(F_mp[it.m], F_mp[it.m], target_prec);

        if(!arb_overlaps(F_mp[it.m], vref_mp))
        {
            std::cout << "Entry failed test: m = " << it.m << " t = " << it.t << "\n";
            char * s1 = arb_get_str(F_mp[it.m], 1000, 0);
            char * s2 = arb_get_str(vref_mp, 1000, 0);
            std::cout << "   Calculated: " << s1 << "\n";
            std::cout << "    Reference: " << s2 << "\n";
            free(s1);
            free(s2);
            nfailed++;
        }
    }

    arb_clear(t_mp);
    arb_clear(vref_mp);
    for(int i = 0; i <= maxm; i++)
        arb_clear(F_mp[i]);

    return nfailed;
}

static long run_boys_test_double(const boys_data & data)
{
    int maxm = 0;
    for(auto & it : data.values)
        maxm = std::max(maxm, it.m);

    std::vector<double> F_dbl(maxm+1);

    long nfailed = 0;

    for(const auto & it : data.values)
    {
        double t_dbl = std::strtod(it.t.c_str(), nullptr);
        double vref_dbl = std::strtod(it.value.c_str(), nullptr);

        mirp_boys_double(F_dbl.data(), it.m, t_dbl);

        if(vref_dbl != F_dbl[it.m])
        {
            std::cout << "Entry failed test: m = " << it.m << " t = " << it.t << "\n";
            std::cout << "  Calculated: " << F_dbl[it.m] << "\n";
            std::cout << "   Reference: " << vref_dbl << "\n";
            nfailed++;
        }
    }

    return nfailed;
}

long run_boys_test(const std::string & filename, const std::string & floattype, long target_prec, long working_prec)
{
    boys_data data;
    try {
        data = read_boys_file(filename);
    }
    catch(std::exception & ex)
    {
        std::cout << "Unable to read data from file \"" << filename << "\"\n";
        std::cout << ex.what() << "\n";
        return 2;
    }

    std::cout << "Read " << data.values.size() << " values from " << filename << "\n";

    long nfailed = 0;

    if(floattype == "mp")
        nfailed = run_boys_test_mp(data, target_prec, working_prec);
    else if(floattype == "interval")
        nfailed = run_boys_test_interval(data, target_prec, working_prec);
    else if(floattype == "double")
        nfailed = run_boys_test_double(data);
    else
    {
        std::string errmsg;
        errmsg = "Unknown floating-point type \"" + floattype + "\"";
        throw std::runtime_error(errmsg);
    }


    double percent_passed = 100.0 - (100.0 * (double)nfailed / (double)data.values.size());
    std::cout << nfailed << " / " << data.values.size() << " failed ("
              << percent_passed << "% passed)\n";

    return nfailed;
}

