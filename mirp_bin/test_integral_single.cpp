/*! \file
 *
 * \brief Functions related to testing four-center single integrals
 */

#include "mirp_bin/testfile_io.hpp"
#include "mirp_bin/test_integral.hpp"
#include "mirp_bin/test_common.hpp"
#include "mirp_bin/callback_helper.hpp"

#include <mirp/pragma.h>
#include <mirp/math.h>

#include <cmath>
#include <iostream>

namespace mirp {

template<int N>
void integral_single_create_test(const std::string & input_filepath,
                                 const std::string & output_filepath,
                                 slong working_prec, long ndigits,
                                 const std::string & header,
                                 typename callback_helper<N>::cb_single_str_type cb)
{
    integral_single_data data = testfile_read_integral_single(input_filepath, N, true);

    data.ndigits = ndigits;
    data.working_prec = working_prec;
    data.header += header;

    /* What we need for the number of digits (plus some safety) */
    const slong min_prec = static_cast<slong>( static_cast<double>(ndigits+5) / MIRP_LOG_10_2 );

    arb_t integral;
    arb_init(integral);

    std::array<std::array<const char *, 3>, N> xyz;
    std::array<std::array<int, 3>, N> lmn;
    std::array<const char *, N> alpha;

    for(auto & ent : data.entries)
    {
        if(ent.g.size() != N)
            throw std::runtime_error("Entry does not have the correct number of gaussians");

        for(int n = 0; n < N; n++)
        {
            lmn[n] = ent.g[n].lmn;
            alpha[n] = ent.g[n].alpha.c_str();

            for(int i = 0; i < 3; i++)
                xyz[n][i] = ent.g[n].xyz[i].c_str();
        }

        callback_helper<N>::call_single_str(integral, lmn, xyz, alpha, working_prec, cb); 

        slong bits = arb_rel_accuracy_bits(integral);
        if(bits > 0 && bits < min_prec)
            throw std::runtime_error("Working precision not large enough for the number of digits");

        char * s = arb_get_str(integral, ndigits, 0);
        ent.integral = s;
        free(s);
    }

    testfile_write_integral_single(output_filepath, data);
    arb_clear(integral);
}

template<int N>
long integral_single_verify_test(const std::string & filepath,
                                 slong working_prec,
                                 typename callback_helper<N>::cb_single_str_type cb)
{
    long nfailed = 0;

    integral_single_data data = testfile_read_integral_single(filepath, N, false);

    arb_t integral, integral_ref;
    arb_init(integral);
    arb_init(integral_ref);

    std::array<std::array<const char *, 3>, N> xyz;
    std::array<std::array<int, 3>, N> lmn;
    std::array<const char *, N> alpha;

    for(const auto & ent : data.entries)
    {
        for(int n = 0; n < N; n++)
        {
            lmn[n] = ent.g[n].lmn;
            alpha[n] = ent.g[n].alpha.c_str();

            for(int i = 0; i < 3; i++)
                xyz[n][i] = ent.g[n].xyz[i].c_str();
        }

        callback_helper<N>::call_single_str(integral, lmn, xyz, alpha, working_prec+16, cb); 

        arb_set_str(integral_ref, ent.integral.c_str(), working_prec);

        /* Rounding the reference value to the working precision results in
         * an interval. Does that interval contain our (more precise) result? */
        if(!arb_overlaps(integral_ref, integral))
        {
            std::cout << "Entry failed test:\n";
            char * s1 = arb_get_str(integral, 2*data.ndigits, 0);
            char * s2 = arb_get_str(integral_ref, 2*data.ndigits, 0);
            std::cout << "   Calculated: " << s1 << "\n";
            std::cout << "    Reference: " << s2 << "\n\n";
            free(s1);
            free(s2);
            nfailed++;
        }
    }

    arb_clear(integral);
    arb_clear(integral_ref);

    print_results(nfailed, data.entries.size());

    return nfailed;
}


template<int N>
long integral_single_verify_test_exact(const std::string & filepath,
                                       typename callback_helper<N>::cb_single_exact_type cb,
                                       typename callback_helper<N>::cb_single_type cb_arb)
{
    long nfailed = 0;

    integral_single_data data = testfile_read_integral_single(filepath, N, false);

    std::array<std::array<int, 3>, N> lmn;
    
    std::array<std::array<double, 3>, N> xyz;
    std::array<double, N> alpha;


    std::array<arb_ptr, N> xyz_arb;
    std::array<arb_t, N> alpha_arb;

    for(auto & it : xyz_arb)
        it = _arb_vec_init(3);
    for(auto & it : alpha_arb)
        arb_init(it);

    arb_t integral_arb;
    arb_init(integral_arb);

    double integral;

    for(const auto & ent : data.entries)
    {
        for(int n = 0; n < N; n++)
        {
            lmn[n] = ent.g[n].lmn;
            alpha[n] = std::strtod(ent.g[n].alpha.c_str(), nullptr);
            arb_set_d(alpha_arb[n], alpha[n]);
            
            for(int i = 0; i < 3; i++)
            {
                xyz[n][i] = std::strtod(ent.g[n].xyz[i].c_str(), nullptr);
                arb_set_d(xyz_arb[n] + i, xyz[n][i]);
            }
        }

        /* compute using the callback */
        callback_helper<N>::call_single_exact(&integral, lmn, xyz, alpha, cb);

        /* Compute using very high precision */
        callback_helper<N>::call_single_arb(integral_arb, lmn, xyz_arb, alpha_arb, 512, cb_arb);

        slong acc_bits = arb_rel_accuracy_bits(integral_arb);

        /* If it's <= 0 that is ok */
        if(acc_bits > 0 && acc_bits < 64)
            throw std::logic_error("Not enough bits in testing exact integral function. Contact the developer");

        double vref_dbl = std::strtod(ent.integral.c_str(), nullptr);
        double vref2_dbl = arf_get_d(arb_midref(integral_arb), ARF_RND_NEAR);

        PRAGMA_WARNING_PUSH
        PRAGMA_WARNING_IGNORE_FP_EQUALITY

        if(integral != vref_dbl && integral != vref2_dbl)
        {
            std::cout << "Entry failed test:\n";
            for(int i = 0; i < N; i++)
            {
                std::cout << ent.g[i].lmn[0] << " "
                          << ent.g[i].lmn[1] << " "
                          << ent.g[i].lmn[2] << " "
                          << ent.g[i].xyz[0] << " "
                          << ent.g[i].xyz[1] << " "
                          << ent.g[i].xyz[2] << " "
                          << ent.g[i].alpha << "\n";
            }

            auto old_cout_prec = std::cout.precision(17);
            std::cout << "     Calculated: " << integral << "\n";
            std::cout << "      Reference: " << vref2_dbl << "\n";
            std::cout << " File Reference: " << vref_dbl << "\n\n";
            std::cout.precision(old_cout_prec);
            nfailed++;
        }

        PRAGMA_WARNING_POP
    }

    arb_clear(integral_arb);
    for(auto & it : xyz_arb)
        _arb_vec_clear(it, 3);
    for(auto & it : alpha_arb)
        arb_clear(it);

    print_results(nfailed, data.entries.size());

    return nfailed;
}


/**********************************
 * Template instantiations
 **********************************/
template void
integral_single_create_test<4>(
        const std::string &, const std::string &,
        slong, long, const std::string &,
        typename callback_helper<4>::cb_single_str_type);

template long
integral_single_verify_test<4>(
        const std::string &, slong,
        callback_helper<4>::cb_single_str_type);

template long
integral_single_verify_test_exact<4>(
        const std::string &,
        callback_helper<4>::cb_single_exact_type,
        callback_helper<4>::cb_single_type);


} // close namespace mirp

