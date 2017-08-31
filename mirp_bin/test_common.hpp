/*! \file
 *
 * \brief Miscellaneous testing functions
 */

#pragma once

#include <string>
#include <vector>
#include <istream>

namespace mirp {


/*! \brief Compare two double-precision numbers for equality within a tolerance
 *
 * This tests the relative difference between two double-precision values.
 *
 * The relative difference is calculated as
 *
 * abs(a - b)/max(abs(a), abs(b))
 *
 * \param [in] a,b Values to compare
 * \param [in] tol Tolerance to compare to
 * \return True if the relative difference is less than \p tol
 */
bool almost_equal(double a, double b, double tol);


/*! \brief Print the results of a test
 *
 * \param [in] nfailed Number of failed tests
 * \param [in] ntests Total number of tests run
 */
void print_results(unsigned long nfailed, unsigned long ntests);


/*! \brief Converts a string to lower case, returning a copy */
std::string to_lower(const std::string & s);


/*! \brief Trims spaces and tabs from both ends of a string, returning a copy */
std::string trim(const std::string & s);


/*! \brief Splits a string into components
 *
 * The components must be separated by a space
 */
std::vector<std::string> split(const std::string & s);


/*! \brief Advances the stream past any comment lines */
void file_skip_comments(std::istream & fs, char commentchar);


} // close namespace mirp

