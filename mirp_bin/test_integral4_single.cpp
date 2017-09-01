/*! \file
 *
 * \brief Functions related to testing four-center single integrals
 */

#include "mirp_bin/testfile_io.hpp"
#include "mirp_bin/test_integral4.hpp"
#include "mirp_bin/test_common.hpp"

#include <mirp/math.h>

#include <cmath>
#include <iostream>

namespace mirp {


void integral4_single_create_test(const std::string & input_filepath,
                                  const std::string & output_filepath,
                                  long ndigits, const std::string & header,
                                  cb_integral4_single_target_str cb)
{
    integral_single_data data = testfile_read_integral_single(input_filepath, 4, true);



    data.ndigits = ndigits;
    data.header += header;

    arb_t integral;
    arb_init(integral);

    const slong target_prec = (ndigits+8) / MIRP_LOG_10_2;

    /* Needed to unpack XYZ */
    const char * A[3];
    const char * B[3];
    const char * C[3];
    const char * D[3];

    for(auto & ent : data.entries)
    {
        if(ent.g.size() != 4)
            throw std::runtime_error("Entry does not have 4 gaussians");

        for(int i = 0; i < 3; i++)
        {
            A[i] = ent.g[0].xyz[i].c_str();
            B[i] = ent.g[1].xyz[i].c_str();
            C[i] = ent.g[2].xyz[i].c_str();
            D[i] = ent.g[3].xyz[i].c_str();
        }

        cb(integral,
           ent.g[0].lmn.data(), A, ent.g[0].alpha.c_str(),
           ent.g[1].lmn.data(), B, ent.g[1].alpha.c_str(),
           ent.g[2].lmn.data(), C, ent.g[2].alpha.c_str(),
           ent.g[3].lmn.data(), D, ent.g[3].alpha.c_str(),
           target_prec);

        char * s = arb_get_str(integral, ndigits, ARB_STR_NO_RADIUS);
        ent.integral = s;
        free(s);
    }

    testfile_write_integral_single(output_filepath, data);
    arb_clear(integral);
}


long integral4_single_run_test(const std::string & filepath,
                               long target_prec,
                               cb_integral4_single_target_str cb)
{
    long nfailed = 0;

    integral_single_data data = testfile_read_integral_single(filepath, 4, false);

    arb_t integral, integral_ref;
    arb_init(integral);
    arb_init(integral_ref);

    /* Number of binary digits contained in the reference value strings
       (and the number of binary digits of accuracy, taking into account
       that the number printed is +/- 1 decimal ulp) */
    const long integral_bits = data.ndigits / MIRP_LOG_10_2;
    const long round_bits = (data.ndigits - 1) / MIRP_LOG_10_2;

    /* Needed to unpack XYZ */
    const char * A[3];
    const char * B[3];
    const char * C[3];
    const char * D[3];

    for(const auto & ent : data.entries)
    {
        for(int i = 0; i < 3; i++)
        {
            A[i] = ent.g[0].xyz[i].c_str();
            B[i] = ent.g[1].xyz[i].c_str();
            C[i] = ent.g[2].xyz[i].c_str();
            D[i] = ent.g[3].xyz[i].c_str();
        }

        cb(integral,
           ent.g[0].lmn.data(), A, ent.g[0].alpha.c_str(),
           ent.g[1].lmn.data(), B, ent.g[1].alpha.c_str(),
           ent.g[2].lmn.data(), C, ent.g[2].alpha.c_str(),
           ent.g[3].lmn.data(), D, ent.g[3].alpha.c_str(),
           target_prec + 16);

        /* Convert the reference to arb_t. The reference is printed to +/- 1 ulp (decimal). */
        if(ent.integral == "0")
            arb_zero(integral_ref);
        else
        {
            arb_set_str(integral_ref, ent.integral.c_str(), integral_bits + 16);
            arf_mag_add_ulp(arb_radref(integral_ref),
                            arb_radref(integral_ref),
                            arb_midref(integral_ref),
                            round_bits);
            arb_set_round(integral_ref, integral_ref, target_prec);
        }

        /* Rounding the reference value to the target precision results in
         * an interval. Does that interval contain our (more precise) result? */
        if(!arb_equal(integral_ref, integral) && !arb_contains(integral_ref, integral))
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


long integral4_single_run_test_d(const std::string & filepath,
                                 cb_integral4_single_d cb)
{
    long nfailed = 0;

    integral_single_data data = testfile_read_integral_single(filepath, 4, false);

    /* Needed to unpack XYZ */
    double A[3];
    double B[3];
    double C[3];
    double D[3];

    for(const auto & ent : data.entries)
    {
        for(int i = 0; i < 3; i++)
        {
            A[i] = std::strtod(ent.g[0].xyz[i].c_str(), nullptr);
            B[i] = std::strtod(ent.g[1].xyz[i].c_str(), nullptr);
            C[i] = std::strtod(ent.g[2].xyz[i].c_str(), nullptr);
            D[i] = std::strtod(ent.g[3].xyz[i].c_str(), nullptr);
        }

        double integral;
        const double alpha1_d = std::strtod(ent.g[0].alpha.c_str(), nullptr);
        const double alpha2_d = std::strtod(ent.g[1].alpha.c_str(), nullptr);
        const double alpha3_d = std::strtod(ent.g[2].alpha.c_str(), nullptr);
        const double alpha4_d = std::strtod(ent.g[3].alpha.c_str(), nullptr);

        cb(&integral,
           ent.g[0].lmn.data(), A, alpha1_d,
           ent.g[1].lmn.data(), B, alpha2_d,
           ent.g[2].lmn.data(), C, alpha3_d,
           ent.g[3].lmn.data(), D, alpha4_d);

        double integral_ref = std::strtod(ent.integral.c_str(), nullptr);

        if(!almost_equal(integral, integral_ref, 1e-13))
        {
            double reldiff = std::fabs(integral_ref - integral);
            reldiff /= std::fmax(std::fabs(integral_ref), std::fabs(integral));

            std::cout << "Entry failed test:\n";
            for(int i = 0; i < 4; i++)
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
            std::cout << "   Calculated: " << integral << "\n";
            std::cout << "    Reference: " << integral_ref << "\n";
            std::cout << "Relative Diff: " << reldiff << "\n\n";
            std::cout.precision(old_cout_prec);
            nfailed++;
        }
    }

    print_results(nfailed, data.entries.size());

    return nfailed;
}


long integral4_single_run_test_exact(const std::string & filepath,
                                     cb_integral4_single_exact cb,
                                     cb_integral4_single_target cb_mp)
{
    long nfailed = 0;

    integral_single_data data = testfile_read_integral_single(filepath, 4, false);

    /* Needed to unpack XYZ */
    double A[3];
    double B[3];
    double C[3];
    double D[3];

    arb_ptr A_mp = _arb_vec_init(3);
    arb_ptr B_mp = _arb_vec_init(3);
    arb_ptr C_mp = _arb_vec_init(3);
    arb_ptr D_mp = _arb_vec_init(3);
    arb_ptr alpha_mp = _arb_vec_init(4);

    arb_t integral_mp;
    arb_init(integral_mp);

    for(const auto & ent : data.entries)
    {
        for(int i = 0; i < 3; i++)
        {
            A[i] = std::strtod(ent.g[0].xyz[i].c_str(), nullptr);
            B[i] = std::strtod(ent.g[1].xyz[i].c_str(), nullptr);
            C[i] = std::strtod(ent.g[2].xyz[i].c_str(), nullptr);
            D[i] = std::strtod(ent.g[3].xyz[i].c_str(), nullptr);
            arb_set_d(A_mp + i, A[i]);
            arb_set_d(B_mp + i, B[i]);
            arb_set_d(C_mp + i, C[i]);
            arb_set_d(D_mp + i, D[i]);
        }

        /* compute using the callback */
        double integral;
        const double alpha1_d = std::strtod(ent.g[0].alpha.c_str(), nullptr);
        const double alpha2_d = std::strtod(ent.g[1].alpha.c_str(), nullptr);
        const double alpha3_d = std::strtod(ent.g[2].alpha.c_str(), nullptr);
        const double alpha4_d = std::strtod(ent.g[3].alpha.c_str(), nullptr);
        arb_set_d(alpha_mp + 0, alpha1_d);
        arb_set_d(alpha_mp + 1, alpha2_d);
        arb_set_d(alpha_mp + 2, alpha3_d);
        arb_set_d(alpha_mp + 3, alpha4_d);

        cb(&integral,
           ent.g[0].lmn.data(), A, alpha1_d,
           ent.g[1].lmn.data(), B, alpha2_d,
           ent.g[2].lmn.data(), C, alpha3_d,
           ent.g[3].lmn.data(), D, alpha4_d);

        /* Compute using very high precision */
        cb_mp(integral_mp,
              ent.g[0].lmn.data(), A_mp, alpha_mp + 0,
              ent.g[1].lmn.data(), B_mp, alpha_mp + 1,
              ent.g[2].lmn.data(), C_mp, alpha_mp + 2,
              ent.g[3].lmn.data(), D_mp, alpha_mp + 3,
              512);


        slong acc_bits = arb_rel_accuracy_bits(integral_mp);

        /* If it's <= 0 that is ok */
        if(acc_bits > 0 && acc_bits < 64)
            throw std::logic_error("Not enough bits in testing exact integral function. Contact the developer");

        double vref_dbl = std::strtod(ent.integral.c_str(), nullptr);
        double vref2_dbl = arf_get_d(arb_midref(integral_mp), ARF_RND_NEAR);

        if(integral != vref_dbl && integral != vref2_dbl)
        {
            std::cout << "Entry failed test:\n";
            for(int i = 0; i < 4; i++)
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
    }

    arb_clear(integral_mp);
    _arb_vec_clear(A_mp, 3);
    _arb_vec_clear(B_mp, 3);
    _arb_vec_clear(C_mp, 3);
    _arb_vec_clear(D_mp, 3);
    _arb_vec_clear(alpha_mp, 4);

    print_results(nfailed, data.entries.size());

    return nfailed;
}

} // close namespace mirp

