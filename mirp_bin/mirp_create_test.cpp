#include "mirp_bin/cmdline.hpp"
#include <iostream>

using namespace mirp;

namespace mirp {

// in boys_create_test.cpp
void boys_create_test(const std::string & infile, const std::string & outfile, long ndigits, const std::string & header);

}


int main(int argc, char ** argv)
{
    std::string infile, outfile;
    std::string integral;
    long ndigits;
    
    try {
        auto cmdline = convert_cmdline(argc, argv);

        infile = cmdline_get_arg_str(cmdline, "--infile");
        outfile = cmdline_get_arg_str(cmdline, "--outfile");
        integral = cmdline_get_arg_str(cmdline, "--integral");
        ndigits = cmdline_get_arg_long(cmdline, "--ndigits");
    }
    catch(std::exception & ex)
    {
        std::cout << "Error parsing command line: " << ex.what() << "\n";
        return 1;
    }

    // Create a header from the command line
    std::string header("# Reference values for the ");
    header += integral;
    header += " integral generated with:\n";
    header += "#  ";
    for(int i = 0; i < argc; i++)
        header += " " + std::string(argv[i]);
    header += "\n#\n";

    try
    {
        if(integral == "boys")
            boys_create_test(infile, outfile, ndigits, header);
        else
        {
            std::cout << "Integral \"" << integral << "\" is not valid\n";
            return 3;
        }
            
        return 0;
    }
    catch(std::exception & ex)
    {
        std::cout << "Error while running tests: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
