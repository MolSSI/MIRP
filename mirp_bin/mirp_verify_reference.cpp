/*! \file
 *
 * \brief mirp_verify_reference main function
 */

#include "mirp_bin/cmdline.hpp"
#include "mirp_bin/ref_integral.hpp"

#include <mirp/kernels/all.h>

#include <sstream>
#include <iostream>
#include <stdexcept>

using namespace mirp;


static void print_help(void)
{
    std::cout << "\n"
              << "mirp_verify_reference - Tests the values in a reference data file\n"
              << "\n"
              << "\n"
              << "Required arguments:\n"
              << "    --file         Reference file to test\n"
              << "    --integral     The type of integral to compute. Possibilities are:\n"
              << "                       eri\n"
              << "\n"
              << "\n"
              << "Other arguments:\n"
              << "    -h, --help     Display this help screen\n"
              << "\n";
}



/*! \brief Main function */
int main(int argc, char ** argv)
{
    std::string infile, integral;

    try {
        auto cmdline = convert_cmdline(argc, argv);
        if(cmdline.size() == 0 || cmdline_get_switch(cmdline, "-h") || cmdline_get_switch(cmdline, "--help"))
        {
            print_help();
            return 0;
        }

        infile = cmdline_get_arg_str(cmdline, "--file");
        integral = cmdline_get_arg_str(cmdline, "--integral");

        if(cmdline.size() != 0)
        {
            std::stringstream ss;
            ss << "Unknown command line arguments:\n";
            for(const auto & it : cmdline)
                ss << "  " << it << "\n";
            throw std::runtime_error(ss.str());
        }
    }
    catch(std::exception & ex)
    {
        std::cout << "\nError parsing command line: " << ex.what() << "\n\n";
        std::cout << "Run \"mirp_verify_reference -h\" for help\n\n";
        return 1;
    }

    try
    {
        long nfailed = -1;

        if(integral == "eri")
        {
            nfailed = integral4_test_reference(infile, mirp_eri_exact);
        }
        else
        {
            std::cout << "Integral \"" << integral << "\" is not valid\n";
            return 3;
        }


        if(nfailed)
            return 1;
        else
            return 0;
    }
    catch(std::exception & ex)
    {
        std::cout << "Error while running tests: " << ex.what() << "\n";
        return 1;
    }
}
