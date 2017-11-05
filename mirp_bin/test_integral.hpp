/*! \file
 *
 * \brief Functions related to testing four-center integrals
 */

#pragma once

#include <mirp/typedefs.h>
#include <string>

#include "mirp_bin/callback_helper.hpp"

namespace mirp {


/************************************************
 * Testing Single Integrals
 ************************************************/

/*! \brief Creates a test of single cartesian integrals
 *
 * Any existing output file (given by \p output_filepath) will be overwritten.
 *
 * \throw std::runtime_error if there is a problem opening the file or there
 *        there is a problem reading or writing the data
 *
 * \throw std::runtime_error if the working precision is not sufficient for
 *        the specified number of digits
 *
 * \tparam N Number of centers the integral needs
 * \param [in] input_filepath  The input file to use for integral parameters
 * \param [in] output_filepath The output file to write the computed integrals to
 * \param [in] working_prec    Internal working precision to use
 * \param [in] ndigits         Number of digits to print
 * \param [in] header          Header information to add to the file
 *                             (appended to the input file header)
 * \param [in] cb              Function that computes single cartesian integrals
 */
template<int N>
void integral_single_create_test(const std::string & input_filepath,
                                 const std::string & output_filepath,
                                 slong working_prec, long ndigits,
                                 const std::string & header,
                                 typename callback_helper<N>::cb_single_str_type cb);


extern template void
integral_single_create_test<4>(
        const std::string &, const std::string &,
        slong, long, const std::string &,
        callback_helper<4>::cb_single_str_type);


/*! \brief Runs a test of single cartesian integrals using interval math
 *
 * \throw std::runtime_error if there is a problem opening the file or there
 *        there is a problem reading or writing the data
 *
 * \tparam N Number of centers the integral needs
 * \param [in] filepath     Path to the file with the reference data
 * \param [in] working_prec Internal working precision to use
 * \param [in] cb           Function that computes single cartesian integrals
 * \return Number of failed tests
 */
template<int N>
long integral_single_verify_test(const std::string & filepath,
                                 slong working_prec,
                                 typename callback_helper<N>::cb_single_str_type cb);

extern template long
integral_single_verify_test<4>(
        const std::string &, slong,
        callback_helper<4>::cb_single_str_type);


/*! \brief Test single cartesian integrals in exact double precision
 *
 * The integrals are tested to be exactly equal to the reference data
 * or to integral computed with very large accuracy.
 *
 * \tparam N Number of centers the integral needs
 * \param [in] filepath  Path to the file with the reference data
 * \param [in] cb        Function that computes single cartesian integrals
 *                       in exact double precision
 * \param [in] cb_arb     Function that computes single cartesian integrals
 *                       using interval arithmetic
 * \return Number of failed tests
 */
template<int N>
long integral_single_verify_test_exact(const std::string & filepath,
                                       typename callback_helper<N>::cb_single_exact_type cb,
                                       typename callback_helper<N>::cb_single_type cb_arb);

extern template long
integral_single_verify_test_exact<4>(
        const std::string &,
        callback_helper<4>::cb_single_exact_type,
        callback_helper<4>::cb_single_type);


/************************************************
 * Testing Contracted Integrals
 ************************************************/

/*! \brief Creates a test of contracted integrals
 *
 * Any existing output file (given by \p output_filepath) will be overwritten.
 *
 * \throw std::runtime_error if there is a problem opening the file or there
 *        there is a problem reading or writing the data
 *
 * \throw std::runtime_error if the working precision is not sufficient for
 *        the specified number of digits
 *
 * \tparam N Number of centers the integral needs
 * \param [in] input_filepath  The input file to use for integral parameters
 * \param [in] output_filepath The output file to write the computed integrals to
 * \param [in] working_prec    Internal working precision to use
 * \param [in] ndigits         Number of digits to print
 * \param [in] header          Header information to add to the file
 *                             (appended to the input file header)
 * \param [in] cb              Function that computes contracted integrals
 */
template<int N>
void integral_create_test(const std::string & input_filepath,
                          const std::string & output_filepath,
                          slong working_prec, long ndigits,
                          const std::string & header,
                          typename callback_helper<N>::cb_str_type cb);

extern template void
integral_create_test<4>(const std::string &,
                        const std::string &,
                        slong, long,
                        const std::string &,
                        callback_helper<4>::cb_str_type);

/*! \brief Runs a test of single cartesian integrals
 *
 * \throw std::runtime_error if there is a problem opening the file or there
 *        there is a problem reading or writing the data
 *
 * \tparam N Number of centers the integral needs
 * \param [in] filepath     Path to the file with the reference data
 * \param [in] working_prec Internal working precision to use
 * \param [in] cb           Function that computes contracted integrals
 * \return Number of failed tests
 */
template<int N>
long integral_verify_test(const std::string & filepath,
                          slong working_prec,
                          typename callback_helper<N>::cb_str_type cb);

extern template long
integral_verify_test<4>(const std::string &, slong,
    callback_helper<4>::cb_str_type);


/*! \brief Test contracted integrals in exact double precision
 *
 * The integrals are tested to be exactly equal to the reference data
 * or to integral computed with very large accuracy.
 *
 * \tparam N Number of centers the integral needs
 * \param [in] filepath  Path to the file with the reference data
 * \param [in] cb        Function that computes contracted integrals
 *                       in exact double precision
 * \param [in] cb_arb     Function that computes contracted integrals
 *                       using interval arithmetic
 * \return Number of failed tests
 */
template<int N>
long integral_verify_test_exact(const std::string & filepath,
                                typename callback_helper<N>::cb_exact_type cb,
                                typename callback_helper<N>::cb_type cb_arb);

extern template long
integral_verify_test_exact<4>(const std::string &,
                              callback_helper<4>::cb_exact_type,
                              callback_helper<4>::cb_type);


} // close namespace mirp

