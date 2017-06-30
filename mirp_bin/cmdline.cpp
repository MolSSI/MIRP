#include "mirp_bin/cmdline.hpp"
#include <algorithm>
#include <stdexcept>
#include <sstream>

namespace mirp {

bool cmdline_has_arg(const std::vector<std::string> & cmdline, const std::string & arg)

{
    auto it = std::find(cmdline.begin(), cmdline.end(), arg);
    return it != cmdline.end();
}

std::string cmdline_get_arg_str(const std::vector<std::string> & cmdline, const std::string & arg)
{
    if(!cmdline_has_arg(cmdline, arg))
    {
        std::string errmsg;
        errmsg = "Required argument " + arg + " is missing";
        throw std::runtime_error(errmsg);
    }

    auto it = std::find(cmdline.begin(), cmdline.end(), arg);
    if(it == cmdline.end())
    {
        std::string errmsg;
        errmsg = "Required argument " + arg + " is missing, but cmdline_has_arg found it. ";
        errmsg += "Please contact the developer";
        throw std::logic_error(errmsg);
    }

    it++;
    if(it == cmdline.end())
    {
        std::string errmsg;
        errmsg = arg + " expects an argument";
        throw std::runtime_error(errmsg);
    }

    return *it; 
}

std::string cmdline_get_arg_str(const std::vector<std::string> & cmdline, const std::string & arg, const std::string & def)
{
    if(!cmdline_has_arg(cmdline, arg))
        return def;
    else
        return cmdline_get_arg_str(cmdline, arg);
}

long cmdline_get_arg_long(const std::vector<std::string> & cmdline, const std::string & arg)
{
    std::string strarg = cmdline_get_arg_str(cmdline, arg);
    std::stringstream ss(strarg); 
    long ret;
    ss >> ret;

    // If the argument is a valid integer, it will be consumed completely
    if(!ss.eof())
    {
        std::string errmsg;
        errmsg = "Argument to " + arg + " is not a valid integer";
        throw std::runtime_error(errmsg);
    }

    return ret;
}

long cmdline_get_arg_long(const std::vector<std::string> & cmdline, const std::string & arg, long def)
{
    if(!cmdline_has_arg(cmdline, arg))
        return def;
    else
        return cmdline_get_arg_long(cmdline, arg);
}

std::vector<std::string> convert_cmdline(int argc, char ** argv)
{
    std::vector<std::string> ret;

    for(int i = 1; i < argc; i++)
    {
        std::string arg(argv[i]);
        size_t equal_pos = arg.find('=');

        if(equal_pos != std::string::npos)
        {
            ret.push_back(arg.substr(0, equal_pos)); 
            ret.push_back(arg.substr(equal_pos+1, std::string::npos));
        }
        else
            ret.push_back(arg);
    }

    return ret;
}

} // close namespace mirp
