/*! \file
 *
 * \brief Helper functions for reading/writing test files
 */

#pragma once

#include <string>

namespace mirp {

// forward declaration
struct integral_single_data;

/*! \brief Read generic single integral test data from a file
 *
 * If \p is_input is set to true, then the returned data
 * does not have the integral_single_data_entry::integral member
 * populated.
 *
 * \throw std::runtime_error if there is a problem opening or
 *        writing to the file
 *
 * \param [in] filepath  Path to the file to read from
 * \param [in] n         Number of centers in the integral (2 center, 4 center, etc)
 * \param [in] is_input  True if the file is a test input file, false if it is a data file
 * \return Data read from the file
 */
integral_single_data integral_single_read_file(const std::string & filepath,
                                               int n, 
                                               bool is_input);


/*! \brief Write generic single integral test data to a file
 *
 * Any existing file at \p filepath will be overwritten
 *
 * \throw std::runtime_error if there is a problem opening or
 *        writing to the file
 *
 * \param [in] filepath  Path to the file to read from
 * \param [in] data      The data to write to the file
 */
void integral_single_write_file(const std::string & filepath,
                                const integral_single_data & data);

} // closing namespace mirp

