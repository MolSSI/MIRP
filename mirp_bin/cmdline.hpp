#pragma once

#include <string>
#include <vector>

namespace mirp {

bool cmdline_has_arg(const std::vector<std::string> & cmdline, const std::string & arg);

std::string cmdline_get_arg_str(const std::vector<std::string> & cmdline, const std::string & arg);
std::string cmdline_get_arg_str(const std::vector<std::string> & cmdline, const std::string & arg, const std::string & def);

long cmdline_get_arg_long(const std::vector<std::string> & cmdline, const std::string & arg);
long cmdline_get_arg_long(const std::vector<std::string> & cmdline, const std::string & arg, long def);

std::vector<std::string> convert_cmdline(int argc, char ** argv);

} // close namespace mirp
