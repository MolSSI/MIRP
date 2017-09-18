/*! \file
 *
 * \brief mirp_create_test main function
 */

#include "mirp_bin/cmdline.hpp"
#include "mirp_bin/test_boys.hpp"
#include "mirp_bin/test_integral4.hpp"

#include <mirp/kernels/all.h>

#include <sstream>
#include <iostream>
#include <stdexcept>

using namespace mirp;


static void print_help(void)
{
    std::cout << "\n"
              << "mirp_create_test - Create a test with reference data for MIRP given an input\n"
              << "\n"
              << "\n"
              << "Required arguments:\n"
              << "    --infile       Input file to use (usually ends in .inp)\n"
              << "    --outfile      Output file (usually ends in .dat). Existing data will\n"
              << "                       be overwritten\n"
              << "    --integral     The type of integral to compute. Possibilities are:\n"
              << "                       boys\n"
              << "                       eri\n"
              << "                       eri_single\n"
              << "    --prec         Working precision to use in the calculation\n"
              << "    --ndigits      Number of decimal digits to write for each integral\n"
              << "\n"
              << "\n"
              << "Other arguments:\n"
              << "    -h, --help     Display this help screen\n"
              << "\n";
}



/*! \brief Main function */
int main(int argc, char ** argv)
{
    std::string infile, outfile;
    std::string integral;
    long ndigits;
    long working_prec;

    try {
        auto cmdline = convert_cmdline(argc, argv);
        if(cmdline.size() == 0 || cmdline_get_switch(cmdline, "-h") || cmdline_get_switch(cmdline, "--help"))
        {
            print_help();
            return 0;
        }

        infile = cmdline_get_arg_str(cmdline, "--infile");
        outfile = cmdline_get_arg_str(cmdline, "--outfile");
        integral = cmdline_get_arg_str(cmdline, "--integral");
        ndigits = cmdline_get_arg_long(cmdline, "--ndigits");
        working_prec = cmdline_get_arg_long(cmdline, "--prec");

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
        std::cout << "Run \"mirp_create_test -h\" for help\n\n";
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
            boys_create_test(infile, outfile, working_prec, ndigits, header);
        else if(integral == "eri")
        {
            integral4_create_test(infile, outfile,
                                          working_prec, ndigits, header,
                                          mirp_eri_str);
        }
        else if(integral == "eri_single")
        {
            integral4_single_create_test(infile, outfile,
                                                 working_prec, ndigits, header,
                                                 mirp_eri_single_str);
        }
        else
        {
            std::cout << "Integral \"" << integral << "\" is not valid\n";
            return 3;
        }
    }
    catch(std::exception & ex)
    {
        std::cout << "Error while creating tests: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
