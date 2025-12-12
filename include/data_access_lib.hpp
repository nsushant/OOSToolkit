#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>
#include <armadillo>

struct DataFrame
{
    std::unordered_map<std::string, std::vector<std::string>> data;
    std::vector<std::string> headers;

    DataFrame(const std::string &filename)
    {

        std::ifstream file(filename);

        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file: " + filename);
        }

        std::string line;
        // Read header line
        if (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string col;

            while (std::getline(ss, col, ','))
            {
                headers.push_back(col);
                data[col] = {};
            }
        }

        // Read data rows
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string cell;
            size_t i = 0;
            while (std::getline(ss, cell, ','))
            {
                if (i < headers.size())
                    data[headers[i]].push_back(cell);
                i++;
            }
        }
    }

    // Pandas-style column access
    std::vector<std::string> &operator[](const std::string &colname)
    {
        if (data.find(colname) == data.end())
        {
            throw std::invalid_argument("Column not found: " + colname);
        }
        return data[colname];
    }

    // Convert a column to numeric type (like df['x'].astype(float))
    arma::vec getNumeric(const std::string &colname) const
    {
        auto it = data.find(colname);
        if (it == data.end())
        {
            throw std::invalid_argument("Column not found: " + colname);
        }

        const auto &str_col = it->second;
        arma::vec v(str_col.size());

        for (size_t i = 0; i < str_col.size(); ++i)
        {
            try
            {
                v(i) = std::stod(str_col[i]);
            }
            catch (...)
            {
                v(i) = arma::datum::nan;
            }
        }

        return v;
    }

    // Optional: print column names
    void showColumns() const
    {
        std::cout << "Columns: ";
        for (const auto &h : headers)
        {
            std::cout << h << " ";
        }
        std::cout << std::endl;
    }
};

arma::uvec find_idxs_of_match(std::vector<std::string> s, std::string value_to_match);