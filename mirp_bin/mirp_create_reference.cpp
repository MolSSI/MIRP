/*! \file
 *
 * \brief mirp_create_reference main function
 */

#include "mirp_bin/cmdline.hpp"
#include "mirp_bin/ref_integral4.hpp"

#include <mirp/kernels/all.h>

#include <sstream>
#include <iostream>

using namespace mirp;


static void print_help(void)
{
    std::cout << "\n"
              << "mirp_create_reference - Create a reference data file for use with other programs\n"
              << "\n"
              << "\n"
              << "Required arguments:\n"
              << "    --basis        Path to a basis set file (usually ends in .bas)\n"
              << "    --geometry     Path to an XYZ file(usually ends in .xyz)\n"
              << "    --outfile      Output file (usually ends in .dat). Existing data will\n"
              << "                       be overwritten\n"
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
    std::string basfile, xyzfile, outfile;
    std::string integral;

    try {
        auto cmdline = convert_cmdline(argc, argv);
        if(cmdline.size() == 0 || cmdline_get_switch(cmdline, "-h") || cmdline_get_switch(cmdline, "--help"))
        {
            print_help();
            return 0;
        }

        basfile = cmdline_get_arg_str(cmdline, "--basis");
        xyzfile = cmdline_get_arg_str(cmdline, "--geometry");
        outfile = cmdline_get_arg_str(cmdline, "--outfile");
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
        std::cout << "Run \"mirp_create_reference -h\" for help\n\n";
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
        if(integral == "eri")
        {
            integral4_create_reference(xyzfile, basfile, outfile, header,
                                          mirp_eri_exact);
        }
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
