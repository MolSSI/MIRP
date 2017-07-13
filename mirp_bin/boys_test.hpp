/*! \file
 *
 * \brief Functions related to testing the Boys function
 */

#pragma once

#include <arb.h>
#include <vector>
#include <string>

namespace mirp {

    /*! \brief An single entry for a Boys function test */ 
    struct boys_data_entry
    {
        int m;              //!< The order of the Boys function
        std::string t;      //!< The value at which the function is evaluated
        std::string value;  //!< The expected (reference) value
    };

    /*! \brief A collection of data for testing the Boys function */
    struct boys_data
    {
        long ndigits;                         //!< The number of 
        std::string header;                   //!< Header or comments about the test
        std::vector<boys_data_entry> values;  //!< Actual data for the tests
    };


    /*! \brief Finds the maximum value of `m` in a Boys test data structure */
    int boys_max_m(const boys_data & data);


    /*! \brief Read a file with reference data for the Boys function
     *
     * \throw std::runtime_error if the file does not exist or there
     *        there is a problem reading/parsing the file
     *
     * \param [in] filepath Path to the file to read
     * \return The data in the file
     */
    boys_data boys_read_file(const std::string & filepath);


    /*! \brief Read an input file for creating a Boys function test
     *
     * This reads in `m` and `t` values for the Boys function, as
     * well as a descriptive header.
     * 
     * \throw std::runtime_error if the file does not exist or there
     *        there is a problem reading/parsing the file
     *
     * \param [in] filepath Path to the file to read
     * \return The data in the file, with the `values` member empty
     */
    boys_data boys_read_input_file(const std::string & filepath);


    /*! \brief Write a file with reference data for the Boys function
     *
     * Any existing file will be overrwitten
     *
     * \throw std::runtime_error if there is a problem opening the file or there
     *        there is a problem writing the data
     *
     * \param [in] filepath Path to the file to write to
     * \param [in] data     The data to write to the file
     */
    void boys_write_file(const std::string & filepath, const boys_data & data);


    /*! \brief Run a test of the Boys function
     *
     * \param [in] filepath    Path to the file to test
     * \param [in] floattype   Type of floating point to test ("double", for example)
     * \param [in] extra_m     Additional `m` values (used to test recurrence relations)
     * \param [in] target_prec Desired final precision
     * \return The number of tests that have failed   
     */
    long boys_run_test(const std::string & filepath, const std::string & floattype, long extra_m, long target_prec);


    /*! \brief Create a test file for the Boys function from a given input file
     *
     * Any existing output file (given by \p output_filepath) will be overwritten.
     *
     * \throw std::runtime_error if there is a problem opening the file or there
     *        there is a problem reading or writing the data
     *
     * \param [in] input_filepath  Path to the input file
     * \param [in] output_filepath File to write the computed data to
     * \param [in] ndigits         Number of digits to compute
     * \param [in] header          Any descriptive header data
     *                             (will be appended to the existing header in the input file)
     */
    void boys_create_test(const std::string & input_filepath, const std::string & output_filepath,
                          long ndigits, const std::string & header);

} // close namespace mirp

