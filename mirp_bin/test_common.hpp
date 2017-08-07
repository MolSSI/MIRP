/*! \file
 *
 * \brief Functions related to testing single ERI
 */

#pragma once

#include "mirp_bin/file_io.hpp"
#include "mirp_bin/data_entry.hpp"

#include <mirp/math.h>
#include <iostream>
#include <stdexcept>

namespace mirp {
namespace detail {


/*! \brief A function that computes single 4-center integrals with interval arithmetic */
typedef void (*cb_single4_interval)(arb_t,
                                    const int *, arb_srcptr, const arb_t,
                                    const int *, arb_srcptr, const arb_t,
                                    const int *, arb_srcptr, const arb_t,
                                    const int *, arb_srcptr, const arb_t,
                                    slong);


/*! \brief Unpacks the arguments for a call to a function that
 *         computes single 4-center integrals
 */
inline void call_callback(arb_t integral, const int * lmn,
                          arb_srcptr xyz, arb_srcptr alpha,
                          slong working_prec, cb_single4_interval cb)
{
    cb(integral,
       lmn+0, xyz+0, alpha+0,
       lmn+3, xyz+3, alpha+1,
       lmn+6, xyz+6, alpha+2,
       lmn+9, xyz+9, alpha+3,
       working_prec);
}
                      

/*! \brief Computes a single integral to a given target precision
 * 
 * This function converts the entry data (strings) to arb_t and
 * runs the callback function.
 *
 * The accuracy/precision of the computed integral will be at
 * least \p target_prec, and will likely be greater.
 *
 * \throw std::runtime_error if there is an incorrect number of
 *        gaussians in the entry
 *
 * \tparam N             Number of centers for the integral
 * \tparam callback_type The type of the function that computes the integrals
 *
 * \param [out] integral     The computed integral
 * \param [in]  ent          Integral test entry data
 * \param [in]  target_prec  Desired precision of the integral
 * \param [in]  cb           Function to use to compute the integrals
 */
template<unsigned int N, typename callback_type>
void integral_single_interval(arb_t integral,
                              const integral_single_data_entry & ent,
                              slong target_prec,
                              callback_type cb)
{
    if(ent.g.size() != N)
    {
        std::string err("Inconsistent number of gaussians.");
        err += "Expected " + std::to_string(N) + ". Got " + std::to_string(ent.g.size());
        throw std::runtime_error(err);
    }
        
    /* Temporaries */
    int lmn[N*3];

    arb_ptr xyz = _arb_vec_init(N*3);
    arb_ptr alpha = _arb_vec_init(N);

    slong working_prec = target_prec;
    bool sufficient_accuracy = false;

    do {
        working_prec += 16;

        for(unsigned int n = 0; n < N; n++)
        {
            arb_set_str(alpha + n, ent.g[n].alpha.c_str(), working_prec);
            arb_set_str(xyz + (n*3+0), ent.g[n].xyz[0].c_str(), working_prec);
            arb_set_str(xyz + (n*3+1), ent.g[n].xyz[1].c_str(), working_prec);
            arb_set_str(xyz + (n*3+2), ent.g[n].xyz[2].c_str(), working_prec);
            lmn[n*3+0] = ent.g[n].lmn[0];
            lmn[n*3+1] = ent.g[n].lmn[1];
            lmn[n*3+2] = ent.g[n].lmn[2];
        }

        /* Run the call back function using the converted values */
        call_callback(integral, lmn, xyz, alpha, working_prec, cb);
        if(arb_rel_accuracy_bits(integral) >= target_prec)
            sufficient_accuracy = true;
            
    } while(!sufficient_accuracy);

    _arb_vec_clear(xyz, N*3);
    _arb_vec_clear(alpha, N);
}


/*! \brief Computes a single integral in double precision
 * 
 * This function converts the entry data (strings) to double and
 * runs the callback function.
 *
 * \throw std::runtime_error if there is an incorrect number of
 *        gaussians in the entry
 *
 * \tparam N             Number of centers for the integral
 * \tparam callback_type The type of the function that computes the integrals
 *
 * \param [out] integral     The computed integral
 * \param [in]  ent          Integral test entry data
 * \param [in]  cb           Function to use to compute the integrals
 */
template<unsigned int N, typename callback_type>
void integral_single_double(double * integral,
                            const integral_single_data_entry & ent,
                            callback_type cb)
{
    if(ent.g.size() != N)
    {
        std::string err("Inconsistent number of gaussians.");
        err += "Expected " + std::to_string(N) + ". Got " + std::to_string(ent.g.size());
        throw std::runtime_error(err);
    }
        
    /* Temporaries */
    int lmn[N*3];
    double xyz[N*3];
    double alpha[N];

    for(unsigned int n = 0; n < N; n++)
    {
        alpha[n] = atof(ent.g[n].alpha.c_str());
        xyz[n*3+0] = atof(ent.g[n].xyz[0].c_str());
        xyz[n*3+1] = atof(ent.g[n].xyz[1].c_str());
        xyz[n*3+2] = atof(ent.g[n].xyz[2].c_str());
        lmn[n*3+0] = ent.g[n].lmn[0];
        lmn[n*3+1] = ent.g[n].lmn[1];
        lmn[n*3+2] = ent.g[n].lmn[2];
    }

    call_callback(integral, lmn, xyz, alpha, cb);
}

} // close namespace detail


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
template<unsigned int N, typename callback_type>
void integral_single_create_test(const std::string & input_filepath,
                                 const std::string & output_filepath,
                                 long ndigits, const std::string & header,
                                 callback_type cb)
{
    integral_single_data data = integral_single_read_file(input_filepath, N, true);
    data.ndigits = ndigits;
    data.header += header;

    arb_t integral;
    arb_init(integral);

    const slong target_prec = (ndigits+4) / MIRP_LOG_10_2; 

    for(auto & ent : data.values)
    {
        detail::integral_single_interval<4>(integral, ent, target_prec, cb);
        char * s = arb_get_str(integral, ndigits, ARB_STR_NO_RADIUS);
        ent.integral = s;
        free(s);
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
template<unsigned int N, typename callback_type>
long integral_single_run_test_interval(const std::string & filepath,
                                       long target_prec,
                                       callback_type cb)
{
    integral_single_data data = integral_single_read_file(filepath, N, false);

    long nfailed = 0;

    arb_t integral, integral_ref;
    arb_init(integral);
    arb_init(integral_ref);

    for(const auto & ent : data.values)
    {
        detail::integral_single_interval<N>(integral, ent, target_prec, cb); 

        /* The computed precision is guaranteed to be at least target_prec,
           but will likely be greater. Therefore, we have to round the
           calculated value and the reference value to the target precision
         */
        arb_set_round(integral, integral, target_prec);
        arb_set_str(integral_ref, ent.integral.c_str(), target_prec);

        if(!arb_overlaps(integral, integral_ref))
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

    return nfailed;

}

} // close namespace mirp

