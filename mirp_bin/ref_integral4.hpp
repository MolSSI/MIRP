/*! \file
 *
 * \brief Functions related to four-center integral reference files
 */

#pragma once

#include <mirp/typedefs.h>
#include <string>

namespace mirp {


/*! \brief Creates a file with exact double reference values of contracted integrals
 *
 * Any existing output file (given by \p output_filepath) will be overwritten.
 *
 * \throw std::runtime_error if there is a problem opening the file or there
 *        there is a problem reading or writing the data
 *
 * \param [in] xyz_filepath    Path to the XYZ file containing the molecule to use
 * \param [in] basis_filepath  Path to a basis set file to use
 * \param [in] output_filepath The output file to write the computed integrals to
 * \param [in] header          Header information to add to the file
 *                             (appended to the input file header)
 * \param [in] cb              Function that computes contracted integrals
 *                             to exact double precision
 */
void integral4_create_reference(const std::string & xyz_filepath,
                                const std::string & basis_filepath,
                                const std::string & output_filepath,
                                const std::string & header,
                                cb_integral4_exact cb);


/*! \brief Tests a reference file for consistency
 *
 * \throw std::runtime_error if there is a problem opening the file or there
 *        there is a problem reading the data
 *
 * \param [in] ref_filepath    Path to the reference file
 * \param [in] cb              Function that computes contracted integrals
 *                             to exact double precision
 * \return Number of failed tests
 */
long integral4_test_reference(const std::string & ref_filepath,
                              cb_integral4_exact cb);


} // close namespace mirp

