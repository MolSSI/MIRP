/*! \file
 *
 * \brief Functions for parsing the command line given to a program
 */

#pragma once

#include <string>
#include <vector>

namespace mirp {

/*! \brief Check to see if the given command line has an argument
 *
 * \param [in] cmdline The command line to check (should be converted already)
 * \param [in] arg     The argument to check for
 * \return True if \p cmdline contains the argument, false otherwise
 */
bool cmdline_has_arg(const std::vector<std::string> & cmdline, const std::string & arg);


/*! \brief Obtain the value of an argument from the command line
 *         as a string
 *
 * \throw std::runtime_error if the argument key or the value is not found
 *
 * \param [in] cmdline The command line to use (should be converted already)
 * \param [in] arg     The argument key to look up
 * \return The value of the argument given on the command line
 */
std::string cmdline_get_arg_str(const std::vector<std::string> & cmdline, const std::string & arg);


/*! \brief Obtain the value of an argument from the command line
 *         as a string, with a default
 *
 * If the argument key is not given on the command line, the default
 * parameter \p def is returned instead.

 * \param [in] cmdline The command line to use (should be converted already)
 * \param [in] arg     The argument key to look up
 * \param [in] def     A default value of the argument to use if \p arg is not
                       given on the command line
 * \return The value of the argument given on the command line, or the value of \p def
 */
std::string cmdline_get_arg_str(const std::vector<std::string> & cmdline, const std::string & arg, const std::string & def);


/*! \brief Obtain the value of an argument from the command line
 *         as a long integer
 *
 * \copydetails cmdline_get_arg_str(const std::vector<std::string> &, const std::string &)
 */
long cmdline_get_arg_long(const std::vector<std::string> & cmdline, const std::string & arg);


/*! \brief Obtain the value of an argument from the command line
 *         as a long integer, with a default
 *
 * \copydetails cmdline_get_arg_str(const std::vector<std::string> &, const std::string &, const std::string &)
 */
long cmdline_get_arg_long(const std::vector<std::string> & cmdline, const std::string & arg, long def);


/*! \brief Convert the command line passed to a program into a vector of strings
 *
 * After conversion, the command line is able to be used in the other `cmdline_` functions
 *
 * \param [in] argc Number of command line arguments
 * \param [in] argv The command line arguments
 * \return The command line, split into a vector of strings and lightly processed.
 */

std::vector<std::string> convert_cmdline(int argc, char ** argv);

} // close namespace mirp
