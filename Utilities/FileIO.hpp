/*
 *
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Slaven Peles <peles2@llnl.gov>.
 * LLNL-CODE-718378.
 * All rights reserved.
 *
 * This file is part of GridKit™. For details, see github.com/LLNL/GridKit
 * Please also read the LICENSE file.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * - Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the disclaimer below.
 * - Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the disclaimer (as noted below) in the
 *   documentation and/or other materials provided with the distribution.
 * - Neither the name of the LLNS/LLNL nor the names of its contributors may
 *   be used to endorse or promote products derived from this software without
 *   specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
 * SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISINGIN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Lawrence Livermore National Laboratory is operated by Lawrence Livermore
 * National Security, LLC, for the U.S. Department of Energy, National
 * Nuclear Security Administration under Contract DE-AC52-07NA27344.
 *
 * This document was prepared as an account of work sponsored by an agency
 * of the United States government. Neither the United States government nor
 * Lawrence Livermore National Security, LLC, nor any of their employees
 * makes any warranty, expressed or implied, or assumes any legal liability
 * or responsibility for the accuracy, completeness, or usefulness of any
 * information, apparatus, product, or process disclosed, or represents that
 * its use would not infringe privately owned rights. Reference herein to
 * any specific commercial product, process, or service by trade name,
 * trademark, manufacturer, or otherwise does not necessarily constitute or
 * imply its endorsement, recommendation, or favoring by the United States
 * government or Lawrence Livermore National Security, LLC. The views and
 * opinions of authors expressed herein do not necessarily state or reflect
 * those of the United States government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product
 * endorsement purposes.
 *
 */

/**
 * @file FileIO.hpp
 * @author Slaven Peles <slaven.peles@pnnl.gov>
 *
 * Contains definition of a utility for reading lookup tables.
 *
 */
#pragma once

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <regex>
#include <sstream>
#include <vector>

#include "MatPowerUtils.hpp"

namespace
{

using namespace GridKit;
using namespace GridKit::MatPowerUtils;

static const std::string matlab_syntax_error{
    "Only a subset of Matlab syntax is supported."
    "\n\t'=' for assignment must be on the same line as the field, eg "
    "`mpc.version = '2'`."
    "\n\tOpen brace ('[') must be on the same line as the field for matrix "
    "initialization."
    "\n\tEach row of a matrix must be terminated by a semicolon."};

std::ostream& logs()
{
#ifndef NDEBUG
  std::cerr << "[FileIO.hpp]: ";
  return std::cerr;
#else
  static std::ofstream ofs;
  ofs.setstate(std::ios_base::badbit);
  return ofs;
#endif
}

void ltrim(std::string& s)
{
  const std::string nothing = "";
  s = std::regex_replace(s, std::regex("^\\s+"), nothing);
}

void rtrim(std::string& s)
{
  const std::string nothing = "";
  s = std::regex_replace(s, std::regex("\\s+$"), nothing);
}

void trim_matlab_comments(std::string& s)
{
  const std::string nothing = "";
  s = std::regex_replace(s, std::regex("%.+"), nothing);
}

// Retrive MATPOWER component from assignment line.
//
// For example, the string "   mpc.bus =  [ ... ] % Some comment" will
// return the value "bus".
std::string getMatPowerComponent(const std::string& line)
{
  logs() << "Getting matpower component from line\n";
  std::regex pat("mpc.([a-zA-Z]+)\\s*=.+");
  std::smatch matches;
  std::string component;
  if (std::regex_match(line, matches, pat))
  {
    component = matches[1].str();
  }
  else
  {
    throw std::runtime_error(matlab_syntax_error + "\nGot line " + line);
  }
  ltrim(component);
  rtrim(component);
  return component;
}

// Ensure that all of the given line has been consumed, and that the only
// remaining non-whitespace character left in the line is a semicolon.
void checkEndOfMatrixRow(std::istream& is)
{
  std::string rest;
  is >> rest;
  ltrim(rest);
  rtrim(rest);
  if (rest != ";")
    throw std::runtime_error(matlab_syntax_error);
}

template <typename IntT = int, typename RealT = double>
void readMatPowerBusRow(const std::string& row, BusRow<IntT, RealT>& br, LoadRow<IntT, RealT>& lr)
{
  logs() << "Parsing MATPOWER bus row\n";
  std::stringstream is(row);
  is >> br.bus_i
     >> br.type
     >> lr.Pd
     >> lr.Qd
     >> br.Gs
     >> br.Bs
     >> br.area
     >> br.Vm
     >> br.Va
     >> br.baseKV
     >> br.zone
     >> br.Vmax
     >> br.Vmin;

  lr.bus_i = br.bus_i;

  // std::cout << br.str();
  // logs() << "Read BusRow with the following values:\n" << br.str();
  // return br;
}

template <typename IntT = int, typename RealT = double>
void readMatPowerGenRow(GenRow<IntT, RealT>& gr, std::string& row)
{
  logs() << "Parsing MATPOWER gen row\n";
  std::stringstream is(row);
  is >> gr.bus >> gr.Pg >> gr.Qg >> gr.Qmax >> gr.Qmin >> gr.Vg >> gr.mBase
     >> gr.status >> gr.Pmax >> gr.Pmin >> gr.Pc1 >> gr.Pc2 >> gr.Qc1min
     >> gr.Qc1max >> gr.Qc2min >> gr.Qc2max >> gr.ramp_agc >> gr.ramp_10
     >> gr.ramp_30 >> gr.ramp_q >> gr.apf;
  checkEndOfMatrixRow(is);
}

template <typename IntT = int, typename RealT = double>
void readMatPowerBranchRow(BranchRow<IntT, RealT>& br, std::string& row)
{
  logs() << "Parsing MATPOWER branch row\n";
  std::stringstream is(row);
  is >> br.fbus >> br.tbus >> br.r >> br.x >> br.b >> br.rateA >> br.rateB
     >> br.rateC >> br.ratio >> br.angle >> br.status >> br.angmin
     >> br.angmax;
  checkEndOfMatrixRow(is);
}

template <typename IntT = int, typename RealT = double>
void readMatPowerGenCostRow(GenCostRow<IntT, RealT>& gcr, std::string& row)
{
  logs() << "Parsing MATPOWER gen cost row\n";
  // Ensure last character is semicolon.
  rtrim(row);
  if (row[row.size() - 1] != ';')
    throw std::runtime_error(matlab_syntax_error + "\nGot line " + row);

  std::stringstream is(row);
  is >> gcr.kind >> gcr.startup >> gcr.shutdown >> gcr.n;

  for (RealT r; is >> r;) {
    gcr.rest.push_back(r);
  }
}

template <typename IntT = int, typename RealT = double>
void readMatPowerVersion(MatPower<IntT, RealT>& mp, std::string& line)
{
  logs() << "Parsing matpower version\n";
  std::regex pat("mpc\\.version\\s*=\\s*'([0-9])';");
  std::smatch matches;
  if (std::regex_match(line, matches, pat)) {
    mp.version = matches[1].str();
  } else {
    throw std::runtime_error(matlab_syntax_error + "\nGot line '" + line + "'");
  }
}

template <typename IntT = int, typename RealT = double>
void readMatPowerBaseMVA(MatPower<IntT, RealT>& mp, std::string& line)
{
  std::regex pat("mpc\\.baseMVA\\s*=\\s*([0-9]+);");
  std::smatch matches;
  if (std::regex_match(line, matches, pat)) {
    std::string s = matches[1];
    mp.baseMVA = std::atoi(s.c_str());
  } else {
    throw std::runtime_error(matlab_syntax_error + "\nGot line '" + line + "'");
  }
}

}  // namespace

namespace GridKit
{
template <typename ScalarT>
void setLookupTable(std::vector<std::vector<ScalarT>>& table,
                    std::string filename,
                    ScalarT& ti,
                    ScalarT& tf)
{
  std::ifstream idata(filename);
  std::string line;
  int oldwordcount = -1;
  while (std::getline(idata, line)) {
    std::istringstream iss(line);
    double word;
    int wordcount = 0;
    std::vector<double> row;
    while (iss >> word) {
      row.push_back(word);
      ++wordcount;
    }
    table.push_back(std::move(row));
    if (oldwordcount != -1) {
      if (oldwordcount != wordcount) {
        std::cerr << "Corrupted input data!\n";
        return;
      }
    } else {
      oldwordcount = wordcount;
    }
  }

  size_t N = table.size();
  ti = table[0][0];
  tf = table[N - 1][0];
}

template <typename ScalarT>
void printLookupTable(std::vector<std::vector<ScalarT>> const& table)
{
  for (size_t i = 0; i < table.size(); ++i) {
    for (size_t j = 0; j < table[i].size(); ++j) {
      std::cout << table[i][j] << " ";
    }
    std::cout << "\n";
  }
}

template <typename IntT = int, typename RealT = double>
void readMatPowerFile(MatPower<IntT, RealT>& mp, std::string& filename)
{
  std::ifstream ifs{filename};
  readMatPower(mp, ifs);
}

template <typename IntT = int,
          typename RealT = double,
          std::size_t MaxLineSize = 1028>
void readMatPower(MatPower<IntT, RealT>& mp, std::istream& is)
{
  using BusRowT     = BusRow<IntT, RealT>;
  using GenRowT     = GenRow<IntT, RealT>;
  using BranchRowT  = BranchRow<IntT, RealT>;
  using GenCostRowT = GenCostRow<IntT, RealT>;
  using LoadRowT    = LoadRow<IntT, RealT>;

  for (std::string line; std::getline(is, line);) {
    // Trim whitespace and remove comments
    ltrim(line);
    rtrim(line);
    logs() << line << "\n";
    trim_matlab_comments(line);

    // Skip empty lines and comment-only lines
    if (line.size() == 0) continue;

    // Skip the matlab function declaration
    if (line.find("function") != std::string::npos) continue;

    // Check for MATPOWER component definitions
    if (line.find("mpc") != std::string::npos) {
      const std::string component = getMatPowerComponent(line);
      logs() << "Got component: '" << component << "'\n";
      // First, parse matrix components
      if (component == "bus") {
        while (std::getline(is, line)) {
          if (line.find("];") != std::string::npos) break;
          BusRowT br;
          LoadRowT lr;
          readMatPowerBusRow(line, br, lr);
          mp.bus.push_back(std::move(br));
          mp.load.push_back(std::move(lr));
        }
      } else if (component == "gen") {
        while (std::getline(is, line)) {
          if (line.find("];") != std::string::npos) break;
          GenRowT gr;
          readMatPowerGenRow(gr, line);
          mp.gen.push_back(gr);
        }
      } else if (component == "branch") {
        while (std::getline(is, line)) {
          if (line.find("];") != std::string::npos) break;
          BranchRowT br;
          readMatPowerBranchRow(br, line);
          mp.branch.push_back(br);
        }
      } else if (component == "gencost") {
        while (std::getline(is, line)) {
          if (line.find("];") != std::string::npos) break;
          GenCostRowT gcr;
          readMatPowerGenCostRow(gcr, line);
          mp.gencost.push_back(gcr);
        }
      }

      // Next, parse scalar components
      else if (component == "version") {
        readMatPowerVersion(mp, line);
      } else if (component == "baseMVA") {
        readMatPowerBaseMVA(mp, line);
      }
    }
  }
}

}  // namespace GridKit
