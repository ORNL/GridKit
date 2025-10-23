
/**
 * @file FileIO.hpp
 * @author Slaven Peles <slaven.peles@pnnl.gov>
 *
 * Contains definition of a utility for reading lookup tables.
 *
 */
#pragma once

#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <vector>

namespace GridKit
{
  /**
   * @brief Reads in an input stream of tabulated data
   *
   * @todo needs to return int for file error codes
   *
   * @tparam ScalarT
   * @param[out] table object in memory where the data from the input stream is
   * @param[in] filename input stream to space and newline separated data
   * @param[out] ti initial time returned
   * @param[out] tf final time returned
   *
   * @pre Input stream should read space separated data. Rows are separated
   * by new line. Each row od data must have the same number of entries. The
   * first column of the data represents time and other columns time dependent
   * variables.
   */
  template <typename ScalarT>
  void setLookupTable(std::vector<std::vector<ScalarT>>& table,
                      std::istream&                      idata,
                      ScalarT&                           ti,
                      ScalarT&                           tf)
  {
    std::string line;
    int         oldwordcount = -1;
    while (std::getline(idata, line))
    {
      std::istringstream  iss(line);
      double              word;
      int                 wordcount = 0;
      std::vector<double> row;
      while (iss >> word)
      {
        row.push_back(word);
        ++wordcount;
      }
      table.push_back(std::move(row));
      if (oldwordcount != -1)
      {
        if (oldwordcount != wordcount)
        {
          std::cerr << "Corrupted input data!\n";
          return;
        }
      }
      else
      {
        oldwordcount = wordcount;
      }
    }

    size_t N = table.size();
    if (N == 0)
    {
      // do nothing
    }
    else if (N == 1)
    {
      ti = table[0][0];
      tf = table[0][0];
    }
    else
    {
      ti = table[0][0];
      tf = table[N - 1][0];
    }
  }

  template <typename ScalarT>
  void printLookupTable(std::vector<std::vector<ScalarT>> const& table)
  {
    for (size_t i = 0; i < table.size(); ++i)
    {
      for (size_t j = 0; j < table[i].size(); ++j)
      {
        std::cout << table[i][j] << " ";
      }
      std::cout << "\n";
    }
  }
} // namespace GridKit
