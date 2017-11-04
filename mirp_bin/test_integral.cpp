/*! \file
 *
 * \brief Functions related to testing four-center integrals
 */

#include "mirp_bin/callback_helper.hpp"
#include "mirp_bin/testfile_io.hpp"
#include "mirp_bin/test_integral.hpp"
#include "mirp_bin/test_common.hpp"

#include <mirp/pragma.h>
#include <mirp/math.h>
#include <mirp/shell.h>

#include <cmath>
#include <iostream>
#include <fstream>

namespace mirp {

namespace detail {

/*! \brief The number of integrals computed in an entry */
static size_t nintegrals(const integral_data_entry & ent)
{
    size_t nint = 1;
    for(const auto & g : ent.g)
        nint *= MIRP_NCART(g.am)*g.ngeneral;
    return nint;
}


template<int N, typename Func>
void integral_create_test(const std::string & input_filepath,
                          const std::string & output_filepath,
                          slong working_prec, long ndigits,
                          const std::string & header,
                          Func cb)
{
    integral_data data = testfile_read_integral(input_filepath, N, true);

    data.ndigits = ndigits;
    data.working_prec = working_prec;
    data.header += header;

    /* What we need for the number of digits (plus some safety) */
    const slong min_prec = static_cast<slong>( static_cast<double>(ndigits+5) / MIRP_LOG_10_2 );

    std::array<std::array<const char *, 3>, N> xyz;
    std::array<std::vector<const char *>, N> alpha, coeff;
    std::array<int, N> am, nprim, ngeneral;

    for(auto & ent : data.entries)
    {
        const size_t nint = nintegrals(ent);
        arb_ptr integrals = _arb_vec_init(nint);

        for(int n = 0; n < N; n++)
        {
            const auto & g = ent.g[n];

            alpha[n].clear();
            coeff[n].clear();

            am[n] = g.am;
            nprim[n] = g.nprim;
            ngeneral[n] = g.ngeneral;

            /* Unpack xyz, exponents, and coefficients */
            for(int i = 0; i < 3; i++)
                xyz[n][i] = g.xyz[i].c_str();
            for(int i = 0; i < g.nprim; i++)
                alpha[n].push_back(g.alpha[i].c_str());
            for(int i = 0; i < g.nprim*g.ngeneral; i++)
                coeff[n].push_back(g.coeff[i].c_str());
        }

        callback_helper<N>::call_str(integrals, am, xyz, nprim, ngeneral, alpha, coeff, working_prec, cb);

        for(size_t i = 0; i < nint; i++)
        {
            slong bits = arb_rel_accuracy_bits(integrals+i);
            if(bits > 0 && bits < min_prec)
                throw std::runtime_error("Working precision not large enough for the number of digits");

            char * s = arb_get_str(integrals+i, ndigits, 0);
            ent.integrals.push_back(s);
            free(s);
        }

        _arb_vec_clear(integrals, nint);
    }

    testfile_write_integral(output_filepath, data);
}


template<int N, typename Func>
long integral_verify_test(const std::string & filepath,
                          slong working_prec,
                          Func cb)
{
    long nfailed = 0;

    integral_data data = testfile_read_integral(filepath, N, false);

    arb_t integral_ref;
    arb_init(integral_ref);

    std::array<std::array<const char *, 3>, N> xyz;
    std::array<std::vector<const char *>, N> alpha, coeff;
    std::array<int, N> am, nprim, ngeneral;

    for(auto & ent : data.entries)
    {
        const size_t nint = nintegrals(ent);
        arb_ptr integrals = _arb_vec_init(nint);

        for(int n = 0; n < N; n++)
        {
            const auto & g = ent.g[n];

            alpha[n].clear();
            coeff[n].clear();

            am[n] = g.am;
            nprim[n] = g.nprim;
            ngeneral[n] = g.ngeneral;

            /* Unpack xyz, exponents, and coefficients */
            for(int i = 0; i < 3; i++)
                xyz[n][i] = g.xyz[i].c_str();
            for(int i = 0; i < g.nprim; i++)
                alpha[n].push_back(g.alpha[i].c_str());
            for(int i = 0; i < g.nprim*g.ngeneral; i++)
                coeff[n].push_back(g.coeff[i].c_str());
        }

        callback_helper<N>::call_str(integrals, am, xyz, nprim, ngeneral, alpha, coeff, working_prec, cb);

        for(size_t i = 0; i < nint; i++)
        {
            arb_set_str(integral_ref, ent.integrals[i].c_str(), working_prec);

            /* Do the intervals overlap? */
            if(!arb_overlaps(integral_ref, integrals+i))
            {
                std::cout << "Entry failed test:\n";
                char * s1 = arb_get_str(integrals+i, 2*data.ndigits, 0);
                char * s2 = arb_get_str(integral_ref, 2*data.ndigits, 0);
                std::cout << "   Calculated: " << s1 << "\n";
                std::cout << "    Reference: " << s2 << "\n\n";
                free(s1);
                free(s2);
                nfailed++;
            }
        }

        _arb_vec_clear(integrals, nint);

    }

    arb_clear(integral_ref);

    print_results(nfailed, data.entries.size());

    return nfailed;
}


template<int N, typename Func, typename Func_arb>
long integral_verify_test_exact(const std::string & filepath,
                                Func cb, Func_arb cb_arb)
{
    long nfailed = 0;

    integral_data data = testfile_read_integral(filepath, N, false);

    std::array<std::array<double, 3>, N> xyz;
    std::array<std::vector<double>, N> alpha, coeff;
    std::vector<double> integrals;

    std::array<arb_ptr, N> xyz_arb;
    for(auto & it : xyz_arb)
        it = _arb_vec_init(3);

    std::array<int, N> am, nprim, ngeneral;

    for(const auto & ent : data.entries)
    {
        const size_t nint = nintegrals(ent);
        integrals.resize(nint);

        arb_ptr integrals_arb = _arb_vec_init(nint);

        std::array<arb_ptr, N> alpha_arb, coeff_arb;

        for(int n = 0; n < N; n++)
        {
            const auto & g = ent.g[n];

            alpha[n].clear();
            coeff[n].clear();

            alpha_arb[n] = _arb_vec_init(g.nprim);
            coeff_arb[n] = _arb_vec_init(g.nprim*g.ngeneral);

            am[n] = g.am;
            nprim[n] = g.nprim;
            ngeneral[n] = g.ngeneral;

            for(int i = 0; i < 3; i++)
            {
                xyz[n][i] = std::strtod(ent.g[n].xyz[i].c_str(), nullptr);
                arb_set_d(xyz_arb[n] + i, xyz[n][i]);
            }

            for(int i = 0; i < g.nprim; i++)
            {
                alpha[n].push_back(std::strtod(g.alpha[i].c_str(), nullptr));
                arb_set_d(alpha_arb[n]+i, alpha[n][i]);
            }

            for(int i = 0; i < g.nprim*g.ngeneral; i++)
            {
                coeff[n].push_back(std::strtod(g.coeff[i].c_str(), nullptr));
                arb_set_d(coeff_arb[n]+i, coeff[n][i]);
            }
        }

        callback_helper<N>::call_exact(integrals.data(), am, xyz, nprim, ngeneral, alpha, coeff, cb);

        /* Compute using very high precision */
        callback_helper<N>::call_arb(integrals_arb, am, xyz_arb, nprim, ngeneral, alpha_arb, coeff_arb, 512, cb_arb);


        for(int n = 0; n < N; n++)
        {
            const auto & g = ent.g[n];
            _arb_vec_clear(alpha_arb[n], g.nprim);
            _arb_vec_clear(coeff_arb[n], g.nprim*g.ngeneral);
        }

        slong acc_bits = mirp_min_accuracy_bits(integrals_arb, nint);

        if(acc_bits > 0 && acc_bits < 64)
            throw std::logic_error("Not enough bits in testing exact integral function. Contact the developer");

        bool failed_shell = false;
        for(size_t i = 0; i < nint; i++)
        {
            double vref_dbl = std::strtod(ent.integrals[i].c_str(), nullptr);
            double vref2_dbl = arf_get_d(arb_midref(integrals_arb+i), ARF_RND_NEAR);

            PRAGMA_WARNING_PUSH
            PRAGMA_WARNING_IGNORE_FP_EQUALITY

            if(integrals[i] != vref_dbl && integrals[i] != vref2_dbl)
            {
                std::cout << "Entry failed test:\n";
                for(int j = 0; j < N; j++)
                {
                    std::cout << ent.g[j].am << " "
                              << ent.g[j].xyz[0] << " "
                              << ent.g[j].xyz[1] << " "
                              << ent.g[j].xyz[2] << "\n";
                }

                auto old_cout_prec = std::cout.precision(17);
                std::cout << "     Calculated: " << integrals[i] << "\n";
                std::cout << "      Reference: " << vref2_dbl << "\n";
                std::cout << " File Reference: " << vref_dbl << "\n\n";
                std::cout.precision(old_cout_prec);
                failed_shell = true;
            }

            PRAGMA_WARNING_POP
        }

        _arb_vec_clear(integrals_arb, nint);

        if(failed_shell)
            nfailed++;
    }

    print_results(nfailed, data.entries.size());

    for(auto & it : xyz_arb)
        _arb_vec_clear(it, 3);

    return nfailed;
}

} // close namespace detail

                        

void integral4_create_test(const std::string & input_filepath,
                           const std::string & output_filepath,
                           slong working_prec, long ndigits,
                           const std::string & header,
                           cb_integral4_str cb)
{
    detail::integral_create_test<4>(input_filepath,
                                    output_filepath,
                                    working_prec, ndigits, header, cb);
}


long integral4_verify_test(const std::string & filepath,
                           slong working_prec,
                           cb_integral4_str cb)
{
    return detail::integral_verify_test<4>(filepath, working_prec, cb);
}

long integral4_verify_test_exact(const std::string & filepath,
                                 cb_integral4_exact cb,
                                 cb_integral4 cb_arb)
{
    return detail::integral_verify_test_exact<4>(filepath, cb, cb_arb);
}




} // close namespace mirp

