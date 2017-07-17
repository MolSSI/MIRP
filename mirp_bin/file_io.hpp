/*! \file
 *
 * \brief Helper functions for reading/writing test files
 */

#pragma once

#include <string>

namespace mirp {

// forward declaration
struct integral_single_data;

integral_single_data integral_single_read_file(const std::string & filepath,
                                               int n, 
                                               bool is_input);


void integral_single_write_file(const std::string & filepath,
                                const integral_single_data & data);

} // closing namespace mirp

