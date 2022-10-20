/*
 *
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Slaven Peles <peles2@llnl.gov>.
 * LLNL-CODE-718378.
 * All rights reserved.
 *
 * This file is part of GridKit. For details, see github.com/LLNL/GridKit
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
 * @file Testing.hpp
 * @author Slaven Peles <slaven.peles@pnnl.gov>
 *
 * Contains utilies for testing.
 *
 */
#pragma once

#include "MatPowerUtils.hpp"
#include <cmath>
#include <iostream>

namespace {

static constexpr double tol_ = 1e-8;

inline std::ostream &errs() {
  std::cerr << "[Utils/Testing.hpp]: ";
  return std::cerr;
}

} // namespace

namespace GridKit {
namespace Testing {

template <typename T> bool isEqual(const T value, const T ref, const T tol) {
  return (std::abs(value - ref) / (1.0 + std::abs(ref)) < tol);
}

template <typename IntT = int, typename RealT = double>
inline bool isEqual(MatPowerUtils::GenCostRow<IntT, RealT> a,
                    MatPowerUtils::GenCostRow<IntT, RealT> b,
                    RealT tol = tol_) {
  int fail = 0;
  fail += a.kind != b.kind;
  fail += a.startup != b.startup;
  fail += a.shutdown != b.shutdown;
  fail += a.n != b.n;
  if (fail) {
    errs() << "Got failure!\na=" << a.str() << "\nb=" << b.str();
  }
  return fail == 0;
}

template <typename IntT = int, typename RealT = double>
inline bool isEqual(MatPowerUtils::GenRow<IntT, RealT> a,
                    MatPowerUtils::GenRow<IntT, RealT> b, RealT tol = tol_) {
  int fail = 0;
  fail += a.bus != b.bus;
  fail += !isEqual(a.Pg, b.Pg, tol);
  fail += !isEqual(a.Qg, b.Qg, tol);
  fail += !isEqual(a.Qmax, b.Qmax, tol);
  fail += !isEqual(a.Qmin, b.Qmin, tol);
  fail += !isEqual(a.Vg, b.Vg, tol);
  fail += a.mBase != b.mBase;
  fail += a.status != b.status;
  fail += a.Pmax != b.Pmax;
  fail += a.Pmin != b.Pmin;
  fail += a.Pc1 != b.Pc1;
  fail += a.Pc2 != b.Pc2;
  fail += a.Qc1min != b.Qc1min;
  fail += a.Qc1max != b.Qc1max;
  fail += a.Qc2min != b.Qc2min;
  fail += a.Qc2max != b.Qc2max;
  fail += a.ramp_agc != b.ramp_agc;
  fail += a.ramp_10 != b.ramp_10;
  fail += a.ramp_30 != b.ramp_30;
  fail += a.ramp_q != b.ramp_q;
  fail += a.apf != b.apf;
  if (fail) {
    errs() << "Got failure!\na=" << a.str() << "\nb=" << b.str();
  }
  return fail == 0;
}

template <typename IntT = int, typename RealT = double>
inline bool isEqual(MatPowerUtils::BusRow<IntT, RealT> a,
                    MatPowerUtils::BusRow<IntT, RealT> b, RealT tol = tol_) {
  int fail = 0;
  fail += a.bus_i != b.bus_i;
  fail += a.type != b.type;
  fail += a.Pd != b.Pd;
  fail += !isEqual(a.Qd, b.Qd, tol);
  fail += a.Gs != b.Gs;
  fail += a.Bs != b.Bs;
  fail += a.area != b.area;
  fail += !isEqual(a.Vm, b.Vm, tol);
  fail += !isEqual(a.Va, b.Va, tol);
  fail += a.baseKV != b.baseKV;
  fail += a.zone != b.zone;
  fail += !isEqual(a.Vmax, b.Vmax, tol);
  fail += !isEqual(a.Vmin, b.Vmin, tol);
  if (fail) {
    errs() << "bus_i: a=" << a.bus_i << ", b=" << b.bus_i << "\n"
           << "type: a=" << a.type << ", b=" << b.type << "\n"
           << "Pd: a=" << a.Pd << ", b=" << b.Pd << "\n"
           << "Qd: a=" << a.Qd << ", b=" << b.Qd << "\n"
           << "Gs: a=" << a.Gs << ", b=" << b.Gs << "\n"
           << "Bs: a=" << a.Bs << ", b=" << b.Bs << "\n"
           << "area: a=" << a.area << ", b=" << b.area << "\n"
           << "Vm: a=" << a.Vm << ", b=" << b.Vm << "\n"
           << "Va: a=" << a.Va << ", b=" << b.Va << "\n"
           << "baseKV: a=" << a.baseKV << ", b=" << b.baseKV << "\n"
           << "zone: a=" << a.zone << ", b=" << b.zone << "\n"
           << "Vmax: a=" << a.Vmax << ", b=" << b.Vmax << "\n"
           << "Vmin: a=" << a.Vmin << ", b=" << b.Vmin << "\n";
  }
  return fail == 0;
}

template <typename IntT = int, typename RealT = double>
inline bool isEqual(MatPowerUtils::LoadRow<IntT, RealT> a,
                    MatPowerUtils::LoadRow<IntT, RealT> b, RealT tol = tol_) {
  int fail = 0;
  fail += a.bus_i != b.bus_i;
  fail += !isEqual(a.Pd, b.Pd, tol);
  fail += !isEqual(a.Qd, b.Qd, tol);
  if (fail) {
    errs() << "bus_i: a=" << a.bus_i << ", b=" << b.bus_i << "\n"
           << "Pd: a=" << a.Pd << ", b=" << b.Pd << "\n"
           << "Qd: a=" << a.Qd << ", b=" << b.Qd << "\n";
  }
  return fail == 0;
}

template <typename IntT = int, typename RealT = double>
inline bool isEqual(MatPowerUtils::BranchRow<IntT, RealT> a,
                    MatPowerUtils::BranchRow<IntT, RealT> b, RealT tol = tol_) {
  int fail = 0;
  fail += a.fbus != b.fbus;
  fail += a.tbus != b.tbus;
  fail += !isEqual(a.r, b.r, tol);
  fail += !isEqual(a.x, b.x, tol);
  fail += !isEqual(a.b, b.b, tol);
  fail += a.rateA != b.rateA;
  fail += a.rateB != b.rateB;
  fail += a.rateC != b.rateC;
  fail += a.ratio != b.ratio;
  fail += a.angle != b.angle;
  fail += a.status != b.status;
  fail += a.angmin != b.angmin;
  fail += a.angmax != b.angmax;
  if (fail) {
    errs() << "Got failure!\na=" << a.str() << "\nb=" << b.str();
  }
  return fail == 0;
}

template <typename T>
inline bool isEqual(std::vector<T> a, std::vector<T> b, double tol = tol_) {
  if (a.size() != b.size())
    throw std::runtime_error([&] {
      std::stringstream errs;
      errs << "Containers do not have the same size!\n"
           << "\tGot a.size() == " << a.size() << "\n"
           << "\tGot b.size() == " << b.size() << "\n";
      return errs.str();
    }());

  int fail = 0;
  for (std::size_t i = 0; i < a.size(); i++) {
    if (!isEqual(a[i], b[i])) {
      fail++;
      errs() << "[isEqual<vector<T>>]: Got failure with i=" << i << ".\n";
    }
  }

  return fail == 0;
}

template <typename IntT = int, typename RealT = double>
inline bool isEqual(MatPowerUtils::MatPower<IntT, RealT> a,
                    MatPowerUtils::MatPower<IntT, RealT> b) {
  int fail = 0;
  fail += a.version != b.version;
  fail += a.baseMVA != b.baseMVA;
  fail += !isEqual(a.bus, b.bus);
  fail += !isEqual(a.gen, b.gen);
  fail += !isEqual(a.gencost, b.gencost);
  fail += !isEqual(a.branch, b.branch);
  return fail == 0;
}

} // namespace Testing

} // namespace GridKit
