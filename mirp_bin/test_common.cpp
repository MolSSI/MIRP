/*! \file
 *
 * \brief Miscellaneous testing functions
 */

#include "test_common.hpp"

#include <mirp/pragma.h>

#include <cmath>
#include <algorithm>
#include <locale> // for std::tolower and std::isspace
#include <map>
#include <iostream>


namespace {

/*! \brief Maps a character to an integer angular momentum ('s' = 0, etc)
 *
 * Angular momentum characters are stored lower case
 */
std::map<char, int> amchar_map
{
    {'s', 0},  {'p', 1},  {'d', 2},  {'f', 3},  {'g', 4},
    {'h', 5},  {'i', 6},  {'k', 7},  {'l', 8},  {'m', 9},
    {'n', 10}, {'o', 11}, {'q', 12}, {'r', 13}, {'t', 14},
    {'u', 15}, {'v', 16}, {'w', 17}, {'x', 18}, {'y', 19},
    {'z', 20}, {'a', 21}, {'b', 22}, {'c', 23}, {'e', 24}
};


/*! \brief Maps an elemental symbol to its atomic Z number ("he" = 2, etc)
 *
 * Element symbols are stored lower case
 */
std::map<std::string, int> z_map
{
    {"h",   1},   {"he", 2},   {"li",  3},   {"be", 4},   {"b",   5},   {"c",   6},   {"n",  7},
    {"o",   8},   {"f", 9},    {"ne",  10},  {"na", 11},  {"mg",  12},  {"al",  13},  {"si", 14},
    {"p",   15},  {"s", 16},   {"cl",  17},  {"ar", 18},  {"k",   19},  {"ca",  20},  {"sc", 21},
    {"ti",  22},  {"v", 23},   {"cr",  24},  {"mn", 25},  {"fe",  26},  {"co",  27},  {"ni", 28},
    {"cu",  29},  {"zn", 30},  {"ga",  31},  {"ge", 32},  {"as",  33},  {"se",  34},  {"br", 35},
    {"kr",  36},  {"rb", 37},  {"sr",  38},  {"y",  39},  {"zr",  40},  {"nb",  41},  {"mo", 42},
    {"tc",  43},  {"ru", 44},  {"rh",  45},  {"pd", 46},  {"ag",  47},  {"cd",  48},  {"in", 49},
    {"sn",  50},  {"sb", 51},  {"te",  52},  {"i",  53},  {"xe",  54},  {"cs",  55},  {"ba", 56},
    {"la",  57},  {"ce", 58},  {"pr",  59},  {"nd", 60},  {"pm",  61},  {"sm",  62},  {"eu", 63},
    {"gd",  64},  {"tb", 65},  {"dy",  66},  {"ho", 67},  {"er",  68},  {"tm",  69},  {"yb", 70},
    {"lu",  71},  {"hf", 72},  {"ta",  73},  {"w",  74},  {"re",  75},  {"os",  76},  {"ir", 77},
    {"pt",  78},  {"au", 79},  {"hg",  80},  {"tl", 81},  {"pb",  82},  {"bi",  83},  {"po", 84},
    {"at",  85},  {"rn", 86},  {"fr",  87},  {"ra", 88},  {"ac",  89},  {"th",  90},  {"pa", 91},
    {"u",   92},  {"np", 93},  {"pu",  94},  {"am", 95},  {"cm",  96},  {"bk",  97},  {"cf", 98},
    {"es",  99},  {"fm", 100}, {"md",  101}, {"no", 102}, {"lr",  103}, {"rf",  104}, {"db", 105},
    {"sg",  106}, {"bh", 107}, {"hs",  108}, {"mt", 109}, {"ds",  110}, {"rg",  111}, {"cn", 112},
    {"uut", 113}, {"fl", 114}, {"uup", 115}, {"lv", 116}, {"uus", 117}, {"uuo", 118}
};

} // close anonymous namespace


namespace mirp {


void print_results(unsigned long nfailed, unsigned long ntests)
{
    double nfailed_d = static_cast<double>(nfailed);
    double ntests_d = static_cast<double>(ntests);
    double percent_passed = 100.0 - (100.0 * nfailed_d / ntests_d);
    std::cout << nfailed << " / " << ntests << " failed ("
              << percent_passed << "% passed)\n";
}


int amchar_to_int(char am)
{
    am = static_cast<char>(std::tolower(am));

    if(amchar_map.find(am) == amchar_map.end())
    {
        std::string err("Unable to convert this character to an angular momentum integer: ");
        err += am;
        throw std::runtime_error(err);
    }
    else
        return amchar_map.at(am);
}


int element_to_z(const std::string & element)
{
    std::string tmp = str_tolower(element);
    if(z_map.find(tmp) == z_map.end())
    {
        std::string err("Unable to convert this element to a Z number: \n");
        err += element;
        throw std::runtime_error(err);
    }
    else
        return z_map.at(tmp);
}


std::string str_tolower(const std::string & s)
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


std::vector<std::string> split(const std::string & s, char sep)
{
    std::vector<std::string> ret;

    const std::string s2 = trim(s);
    const char * end = s2.c_str();

    do
    {
        const char * begin = end;
        while(*begin == sep)
            begin++;

        end = begin;
        while(*end && *end != sep)
            end++;

        ret.push_back(std::string(begin, end));
    }while(0 != *end++);

    return ret;
}


bool file_skip(std::istream & fs, char commentchar)
{
    while(true)
    {
        int c = fs.peek();

        if(!fs.good())                    // did we reach eof?
            break;
        else if(c == commentchar)
            fs.ignore(max_length, '\n');  // skip the rest of the line
        else if(std::isspace(c))          // is whitespace? (includes \n, \t, etc)
            fs.get();                     // discard it
        else
            break;                        // We found some data
    }

    return fs.good();
}


} // close namespace mirp
