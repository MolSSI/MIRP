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
    long ndigits;                         //!< The number of digits of accuracy in the file
    long working_prec;                    //!< Working precision used for the integrals
    std::string header;                   //!< Header or comments about the test
    std::vector<boys_data_entry> entries; //!< Actual data for the tests
};


/*! \brief Finds the maximum value of `m` in a Boys test data structure */
int boys_max_m(const boys_data & data);


/*! \brief Read a file with reference data for the Boys function
 *
 * This reads in `m` and `t` entries for the Boys function, as
 * well as a descriptive header.
 *
 * If \p is_input is true, then then `ndigits` and the reference
 * entries are also read.
 *
 * \throw std::runtime_error if the file does not exist or there
 *        there is a problem reading/parsing the file
 *
 * \param [in] filepath Path to the file to read
 * \param [in] is_input Set to true when reading an input file
 * \return The data in the file
 */
boys_data boys_read_file(const std::string & filepath, bool is_input);


/*! \brief Write a file with reference data for the Boys function
 *
 * This writes reference data files only. It can not be used to write
 * input files.
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
 * \param [in] filepath     Path to the file to test
 * \param [in] floattype    Type of floating point to test ("double", for example)
 * \param [in] extra_m      Additional `m` entries (used to test recurrence relations)
 * \param [in] working_prec Internal working precision to use
 * \return The number of tests that have failed
 */
long boys_verify_test_main(const std::string & filepath,
                           const std::string & floattype,
                           int extra_m, slong working_prec);


/*! \brief Create a test file for the Boys function from a given input file
 *
 * Any existing output file (given by \p output_filepath) will be overwritten.
 *
 * \throw std::runtime_error if there is a problem opening the file or there
 *        there is a problem reading or writing the data
 *
 * \throw std::runtime_error if the working precision is not sufficient for
 *        the specified number of digits
 *
 * \param [in] input_filepath  Path to the input file
 * \param [in] output_filepath File to write the computed data to
 * \param [in] working_prec    Internal working precision to use
 * \param [in] ndigits         Number of decimal digits to compute
 * \param [in] header          Any descriptive header data
 *                             (will be appended to the existing header in the input file)
 */
void boys_create_test(const std::string & input_filepath,
                      const std::string & output_filepath,
                      slong working_prec, long ndigits,
                      const std::string & header);

} // close namespace mirp

