/*! \file
 *
 * \brief Functions related to testing single ERI
 */

#pragma once

#include "mirp_bin/shell.hpp"
#include <arb.h>

namespace mirp {

    struct eri_single_data_entry
    {
        std::array<gaussian_single, 4> g;    //!< Shell for which the data is calculated
        std::string integral; //!< Actual data (integral) for this gaussian
    };

    struct eri_single_data
    {
        long ndigits;                        //!< The number of digits of accuracy in the file 
        std::string header;                  //!< Header or comments about the test
        std::vector<eri_single_data_entry> values;  //!< Actual data for the integrals
    };


    /*! \brief Read a file with reference data for ERI
     *
     * \throw std::runtime_error if the file does not exist or there
     *        there is a problem reading/parsing the file
     *
     * \param [in] filepath Path to the file to read
     * \return The data in the file
     */
    eri_single_data eri_single_read_file(const std::string & filepath);


    /*! \brief Read an input file for creating an ERI test
     *
     * This reads in all the shell information for the ERI, as
     * well as a descriptive header.
     * 
     * \throw std::runtime_error if the file does not exist or there
     *        there is a problem reading/parsing the file
     *
     * \param [in] filepath Path to the file to read
     * \return The data in the file, with the `values` member empty
     */
    eri_single_data eri_single_read_input_file(const std::string & filepath);


    /*! \brief Write a file with reference data for ERI
     *
     * Any existing file will be overrwitten
     *
     * \throw std::runtime_error if there is a problem opening the file or there
     *        there is a problem writing the data
     *
     * \param [in] filepath Path to the file to write to
     * \param [in] data     The data to write to the file
     */
    void eri_single_write_file(const std::string & filepath, const eri_single_data & data);


    /*! \brief Run a test of the ERI functions
     *
     * \param [in] filepath    Path to the file to test
     * \param [in] floattype   Type of floating point to test ("double", for example)
     * \param [in] extra_m     Additional `m` values (used to test recurrence relations)
     * \param [in] target_prec Desired final binary precision (in bits)
     * \return The number of tests that have failed   
     */
    long eri_single_run_test(const std::string & filepath, const std::string & floattype,
                             long extra_m, long target_prec);


    /*! \brief Create a test file for ERI from a given input file
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
     */
    void eri_single_create_test(const std::string & input_filepath,
                                const std::string & output_filepath,
                                long ndigits, const std::string & header);

} // close namespace mirp

