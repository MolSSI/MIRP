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
        std::string header;
        std::vector<boys_data_entry> values;
    };

    boys_data boys_read_file(const std::string & filename);

    boys_data boys_read_input_file(const std::string & filename);

    void boys_write_file(const std::string & filename, const boys_data & data);

    int boys_max_m(const boys_data & data);

} // close namespace mirp

