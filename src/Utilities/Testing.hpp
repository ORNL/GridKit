/**
 * @file Testing.hpp
 * @author Slaven Peles <slaven.peles@pnnl.gov>
 *
 * Contains utilies for testing.
 *
 */
#pragma once

#include <cmath>
#include <iostream>
#include <limits>

#include <PowerSystemData.hpp>

namespace
{

  static constexpr double tol_ = 1e-8;

  inline std::ostream& errs()
  {
    std::cerr << "[Utils/Testing.hpp]: ";
    return std::cerr;
  }

} // namespace

namespace GridKit
{
  namespace Testing
  {

    template <typename T>
    bool isEqual(const T value,
                 const T ref,
                 const T tol = std::numeric_limits<T>::epsilon())
    {
      T error = std::abs(value - ref) / (1.0 + std::abs(ref));
      return (error < tol);
    }

    template <typename RealT = double, typename IdxT = int>
    inline bool isEqual(PowerSystemData::GenCostData<RealT, IdxT> a,
                        PowerSystemData::GenCostData<RealT, IdxT> b,
                        RealT                                     tol = tol_)
    {
      (void) tol; // suppress warning
      int fail  = 0;
      fail     += a.kind != b.kind;
      fail     += a.startup != b.startup;
      fail     += a.shutdown != b.shutdown;
      fail     += a.n != b.n;
      if (fail)
      {
        errs() << "Got failure!\na=" << a.str() << "\nb=" << b.str();
      }
      return fail == 0;
    }

    template <typename RealT = double, typename IdxT = int>
    inline bool isEqual(PowerSystemData::GenData<RealT, IdxT> a,
                        PowerSystemData::GenData<RealT, IdxT> b,
                        RealT                                 tol = tol_)
    {
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

      if (fail)
      {
        errs() << "Got failure!\na=" << a.str() << "\nb=" << b.str();
      }
      return fail == 0;
    }

    template <typename RealT = double, typename IdxT = int>
    inline bool isEqual(PowerSystemData::BusData<RealT, IdxT> a,
                        PowerSystemData::BusData<RealT, IdxT> b,
                        RealT                                 tol = tol_)
    {
      int fail = 0;

      fail += a.bus_i != b.bus_i;
      fail += a.type != b.type;
      fail += a.Gs != b.Gs;
      fail += a.Bs != b.Bs;
      fail += a.area != b.area;
      fail += !isEqual(a.Vm, b.Vm, tol);
      fail += !isEqual(a.Va, b.Va, tol);
      fail += a.baseKV != b.baseKV;
      fail += a.zone != b.zone;
      fail += !isEqual(a.Vmax, b.Vmax, tol);
      fail += !isEqual(a.Vmin, b.Vmin, tol);

      if (fail)
      {
        errs() << "bus_i: a=" << a.bus_i << ", b=" << b.bus_i << "\n"
               << "type: a=" << a.type << ", b=" << b.type << "\n"
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

    template <typename RealT = double, typename IdxT = int>
    inline bool isEqual(PowerSystemData::LoadData<RealT, IdxT> a,
                        PowerSystemData::LoadData<RealT, IdxT> b,
                        RealT                                  tol = tol_)
    {
      int fail = 0;

      fail += a.bus_i != b.bus_i;
      fail += !isEqual(a.Pd, b.Pd, tol);
      fail += !isEqual(a.Qd, b.Qd, tol);

      if (fail)
      {
        errs() << "bus_i: a=" << a.bus_i << ", b=" << b.bus_i << "\n"
               << "Pd: a=" << a.Pd << ", b=" << b.Pd << "\n"
               << "Qd: a=" << a.Qd << ", b=" << b.Qd << "\n";
      }
      return fail == 0;
    }

    template <typename RealT = double, typename IdxT = int>
    inline bool isEqual(PowerSystemData::BranchData<RealT, IdxT> a,
                        PowerSystemData::BranchData<RealT, IdxT> b,
                        RealT                                    tol = tol_)
    {
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

      if (fail)
      {
        errs() << "Got failure!\na=" << a.str() << "\nb=" << b.str();
      }
      return fail == 0;
    }

    template <typename T>
    inline bool isEqual(std::vector<T> a, std::vector<T> b, double tol = tol_)
    {
      if (a.size() != b.size())
        throw std::runtime_error([&]
                                 {
      std::stringstream errs;
      errs << "Containers do not have the same size!\n"
           << "\tGot a.size() == " << a.size() << "\n"
           << "\tGot b.size() == " << b.size() << "\n";
      return errs.str(); }());

      int fail = 0;
      for (std::size_t i = 0; i < a.size(); i++)
      {
        if (!isEqual(a[i], b[i], tol))
        {
          fail++;
          errs() << "[isEqual<vector<T>>]: Got failure with i=" << i << ".\n";
        }
      }

      return fail == 0;
    }

    template <typename RealT = double, typename IdxT = int>
    inline bool isEqual(PowerSystemData::SystemModelData<RealT, IdxT> a,
                        PowerSystemData::SystemModelData<RealT, IdxT> b)
    {
      int fail = 0;

      fail += a.version != b.version;
      fail += a.baseMVA != b.baseMVA;
      fail += !isEqual(a.bus, b.bus);
      fail += !isEqual(a.gen, b.gen);
      fail += !isEqual(a.gencost, b.gencost);
      fail += !isEqual(a.branch, b.branch);
      fail += !isEqual(a.load, b.load);

      return fail == 0;
    }

    template <typename IdxT = size_t, typename RealT = double>
    inline bool isEqual(std::map<IdxT, RealT> a,
                        std::map<IdxT, RealT> b,
                        const RealT tol = std::numeric_limits<RealT>::epsilon())
    {
      int fail = 0;

      fail += a.size() != b.size();

      for (const auto& pair_a : a)
      {
        auto it = b.find(pair_a.first);
        std::cout << pair_a.second << ", " << it->second << std::endl;
        fail += !isEqual(pair_a.second, it->second);
      }

      return fail == 0;
    }

  } // namespace Testing

} // namespace GridKit
