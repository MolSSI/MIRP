/*! \file
 *
 * \brief Functions related to testing single ERI
 */

#include "mirp_bin/file_io.hpp"
#include "mirp_bin/data_entry.hpp"
#include "mirp/typedefs.h"
#include "mirp/math.h"

#include "test_common.hpp"

#include <iostream>

namespace mirp {
namespace detail {


/*! \brief Create a test file for an integral from a given input file
 *
 * Any existing output file (given by \p output_filepath) will be overwritten.
 *
 * \throw std::runtime_error if there is a problem opening the file or there
 *        there is a problem reading or writing the data
 *
 * \param [in] input_filepath  Path to the input file
 * \param [in] output_filepath File to write the computed data to
 * \param [in] ndigits         Number of decimal digits to compute
 * \param [in] header          Any descriptive header data
 *                             (will be appended to the existing header in the input file)
 * \param [in] cb              Function that computes a single integral
 */
void integral4_single_create_test(const std::string & input_filepath,
                                  const std::string & output_filepath,
                                  long ndigits, const std::string & header,
                                  cb_integral4_single_target_str cb)
{
    integral_single_data data = integral_single_read_file(input_filepath, 4, true);

    if(data.values.size() != 4)
    {
        throw std::runtime_error("TODO");
    }


    data.ndigits = ndigits;
    data.header += header;

    arb_t integral;
    arb_init(integral);

    const slong target_prec = (ndigits+4) / MIRP_LOG_10_2; 

    /* We need to unpack XYZ */
    const char * A[3];
    const char * B[3];
    const char * C[3];
    const char * D[3];

    for(auto & ent : data.values)
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
           target_prec);

        
        if(mirp_test_zero_prec(integral, target_prec))
            ent.integral = "0";
        else
        {
            char * s = arb_get_str(integral, ndigits, ARB_STR_NO_RADIUS);
            ent.integral = s;
            free(s);
        }
    }

    integral_single_write_file(output_filepath, data);
    arb_clear(integral);
}


/*! \brief Run tests located in a test file using interval arithmetic
 *
 * \throw std::runtime_error if there is a problem opening the file or there
 *        there is a problem reading or writing the data
 *
 * \tparam N             Number of centers for the integral
 * \tparam callback_type The type of the function that computes the integrals
 *
 * \param [in] filepath    Path to the test data file
 * \param [in] target_prec The target precision (binary precision) to calculate
 * \param [in] cb          Function that computes single integrals
 * \return The number of tests that have failed
 */
long integral4_single_run_test(const std::string & filepath,
                               long target_prec,
                               cb_integral4_single_target_str cb)
{
    long nfailed = 0;

    integral_single_data data = integral_single_read_file(filepath, 4, false);

    arb_t integral, integral_ref;
    arb_init(integral);
    arb_init(integral_ref);

    /* We need to unpack XYZ */
    const char * A[3];
    const char * B[3];
    const char * C[3];
    const char * D[3];

    for(const auto & ent : data.values)
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

        /*
           The computed precision is guaranteed to be at least target_prec,
           but will likely be greater. Round the reference value,
           introducing error. Is the more precise calculated value
           within those error bounds?
         */
        arb_set_str(integral_ref, ent.integral.c_str(), target_prec);

        /* 1.) Test if the calculated value is within the error of the reference
         * 2.) If it's not, test if the calculated value is [0 +/ value] and is
         *     that zero within the target precision
         * 3.) If (2) is true, and the reference integral is exactly zero, then
         *     they are considered equal (and this block is not entered)
         */
        if(!arb_contains(integral_ref, integral) &&
           !(arb_is_zero(integral_ref) && mirp_test_zero_prec(integral, target_prec)))
        {
            std::cout << "Entry failed test:\n";
            char * s1 = arb_get_str(integral, 1000, 0);
            char * s2 = arb_get_str(integral_ref, 1000, 0);
            std::cout << "   Calculated: " << s1 << "\n";
            std::cout << "    Reference: " << s2 << "\n\n";
            free(s1);
            free(s2);
            nfailed++;
        }
    }

    arb_clear(integral);
    arb_clear(integral_ref);

    print_results(nfailed, data.values.size());

    return nfailed;
}


/*! \brief Run tests located in a test file using double precision
 *
 * \throw std::runtime_error if there is a problem opening the file or there
 *        there is a problem reading or writing the data
 *
 * \tparam N             Number of centers for the integral
 * \tparam callback_type The type of the function that computes the integrals
 *
 * \param [in] filepath    Path to the test data file
 * \param [in] target_prec The target precision (binary precision) to calculate
 * \param [in] cb          Function that computes single integrals
 * \return The number of tests that have failed
 */
long integral4_single_run_test_d(const std::string & filepath,
                                 cb_integral4_single_d cb)
{
    long nfailed = 0;

    integral_single_data data = integral_single_read_file(filepath, 4, false);

    /* We need to unpack XYZ */
    double A[3];
    double B[3];
    double C[3];
    double D[3];

    for(const auto & ent : data.values)
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

    print_results(nfailed, data.values.size());

    return nfailed;
}


long integral4_single_run_test_exact(const std::string & filepath,
                                     cb_integral4_single_exact cb,
                                     cb_integral4_single_target cb_mp)
{
    long nfailed = 0;

    integral_single_data data = integral_single_read_file(filepath, 4, false);

    /* We need to unpack XYZ */
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

    for(const auto & ent : data.values)
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

    print_results(nfailed, data.values.size());

    return nfailed;
}

} // close namespace detail
} // close namespace mirp

