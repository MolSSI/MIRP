#pragma once

#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <stdexcept>

namespace mirp {

    struct boys_data_entry
    {
        int m;
        std::string t;
        std::string value;
    };

    struct boys_data
    {
        long ndigits;
        std::vector<boys_data_entry> values;
    };

    boys_data read_boys_file(const std::string & filename);

    void write_boys_file(const std::string & filename, const boys_data & data);

} // close namespace mirp

