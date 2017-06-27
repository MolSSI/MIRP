#include "cmdline.hpp"
#include <iostream>

using namespace mirp;

namespace mirp {

// in boys_run_test.cpp
long boys_run_test(const std::string & filename, const std::string & floattype, long extra_m, long target_prec, long working_prec);

}


int main(int argc, char ** argv)
{
    std::string filename;
    std::string integral;
    std::string floattype;
    long target_prec, working_prec, extra_m;
    
    try {
        auto cmdline = convert_cmdline(argc, argv);

        filename = cmdline_get_arg_str(cmdline, "--filename");
        integral = cmdline_get_arg_str(cmdline, "--integral");
        floattype = cmdline_get_arg_str(cmdline, "--float");

        if(floattype != "double")
        {
            target_prec = cmdline_get_arg_long(cmdline, "--prec");
            working_prec = cmdline_get_arg_long(cmdline, "--working-prec", target_prec*2);
        }

        if(integral == "boys")
            extra_m = cmdline_get_arg_long(cmdline, "--extra-m", 0);
        
    }
    catch(std::exception & ex)
    {
        std::cout << "Error parsing command line: " << ex.what() << "\n";
        return 1;
    }


    try
    {
        long nfailed;
        if(integral == "boys")
            nfailed = boys_run_test(filename, floattype, extra_m, target_prec, working_prec);
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

    return 0;
}
