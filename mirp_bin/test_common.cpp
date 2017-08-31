/*! \file
 *
 * \brief Miscellaneous testing functions
 */

#include "test_common.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>

namespace mirp {

bool almost_equal(double a, double b, double tol)
{
    using std::fabs;
    using std::fmax;

    double diff = fabs(a-b);

    if(diff == 0)
        return true;

    double denominator = fmax(fabs(a), fabs(b));
    return (diff / denominator) < tol;
}


void print_results(unsigned long nfailed, unsigned long ntests)
{
    double percent_passed = 100.0 - (100.0 * (double)nfailed / (double)ntests);
    std::cout << nfailed << " / " << ntests << " failed ("
              << percent_passed << "% passed)\n";
}

std::string to_lower(const std::string & s)
{
    std::string s_copy = s;
    std::transform(s_copy.begin(), s_copy.end(),
                   s_copy.begin(), ::tolower);
    return s_copy;
}


std::string trim(const std::string & s)
{
    const char * ws = " \t";
    const auto start = s.find_first_not_of(ws);
    if(start == std::string::npos)
        return "";

    const auto end = s.find_last_not_of(ws);
    const auto length = end - start + 1;

    return s.substr(start, length);
}


std::vector<std::string> split(const std::string & s)
{
    std::vector<std::string> ret;

    const std::string s2 = trim(s);
    const char * end = s2.c_str();

    do
    {
        const char * begin = end;
        while(*begin == ' ')
            begin++;

        end = begin;
        while(*end && *end != ' ')
            end++;

        ret.push_back(std::string(begin, end));
    }while(0 != *end++);

    return ret;
}


void file_skip_comments(std::istream & fs, char commentchar)
{
    std::string line;
    while(fs.good())
    {
        if(fs.peek() == commentchar)
            std::getline(fs, line);
        else
            break;
    }
}

} // close namespace mirp
