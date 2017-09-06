/*! \file
 *
 * \brief Helper functions for reading/writing reference data files
 */

#pragma once

#include "mirp_bin/data_entry.hpp"
#include <iostream>

namespace mirp {


/*! \brief Reads a double precision number from an ASCII file
 *
 * The main purpose of this is to convert hexfloat to doubles,
 * since it seems streams can't read hexfloat
 *
 * One double is extracted, and the stream advanced
 */
double read_double(std::istream & fs);


/*! \brief Write basis information to a file
 *
 * \param [in] shells The basis information to write
 * \param [in] fs Stream to write the information to
 */
void reffile_write_basis(const std::vector<gaussian_shell> & shells, std::ostream & fs);


/*! \brief Reads basis information from a reference file
 *
 * \param [in] fs The file stream to read from
 * \return The shells contained in the basis in the file
 */
std::vector<gaussian_shell> reffile_read_basis(std::istream & fs);


} // closing namespace mirp

