#pragma once
#include <regex>
#include <string>
#include <iomanip>
#include <sstream>

/**
 *
 * @file PowerSystemData.hpp
 * @author Asher Mancinelli <asher.mancinelli@pnnl.gov>
 *
 * @remark `std::stringstream` is preferred over `operator+(std::string, ...)`
 * since stringstream does not reallocate on append.
 *
 */

namespace GridKit
{

namespace PowerSystemData
{

  template <typename RealT = double, typename IdxT = int>
  struct BusData 
  {
    IdxT  bus_i;  ///< Bus ID
    IdxT  type;   ///< Bus type: 1 = PQ, 2 = PV, 3 = ref, 4 = isolated
    RealT Gs;     ///< Shunt conductance (MW demanded at V = 1.0 p.u.)
    RealT Bs;     ///< Shunt susceptance (MVAr injected at V = 1.0 p.u.)
    IdxT  area;   ///< Area number (>0)
    RealT Vm;     ///< Voltage magnitude (p.u.)
    RealT Va;     ///< Voltage phase (deg)
    RealT baseKV; ///< Base voltage [kV]
    IdxT  zone;   ///< Loss zone number (>0)
    RealT Vmax;   ///< Maximum voltage magnitude (p.u.)
    RealT Vmin;   ///< Minimum voltage magnitude (p.u.)

    inline std::string str() const
    {
      std::stringstream ss;
      std::cerr << std::setw(10) << bus_i
                << std::setw(10) << type 
                << std::setw(10) << Gs
                << std::setw(10) << Bs
                << std::setw(10) << area
                << std::setw(10) << Vm
                << std::setw(10) << Va
                << std::setw(10) << baseKV
                << std::setw(10) << zone
                << std::setw(10) << Vmax
                << std::setw(10) << Vmin;
      ss << "\n";
      return ss.str();
    }
  };

  template <typename RealT = double, typename IdxT = int>
  struct LoadData 
  {
    IdxT  bus_i;  ///< Bus ID
    RealT Pd;     ///< Active power demand [MW]
    RealT Qd;     ///< Reactive power demand [MVAr]

    inline std::string str() const
    {
      std::stringstream ss;
      std::cerr << std::setw(10) << bus_i
                << std::setw(10) << Pd 
                << std::setw(10) << Qd;
      ss << "\n";
      return ss.str();
    }
  };

  template <typename RealT = double, typename IdxT = int>
  struct GenData
  {
    IdxT  bus;      ///< Bus ID
    RealT Pg;       ///< Active power output [MW]
    RealT Qg;       ///< Reactive power output [MVAr]
    RealT Qmax;     ///< Maximum reactive power output [MVAr]
    RealT Qmin;     ///< Minimum reactive power output [MVAr]
    RealT Vg;       ///<
    RealT mBase;    ///< Total MVA base of machine
    IdxT  status;   ///< Service status (>0 in service, <=0 out of service) 
    RealT Pmax;     ///< Maximum active power output [MVAr]
    RealT Pmin;     ///< Minimum active power output [MVAr]
    RealT Pc1;      ///<
    RealT Pc2;      ///<
    RealT Qc1min;   ///<
    RealT Qc1max;   ///<
    RealT Qc2min;   ///<
    RealT Qc2max;   ///<
    RealT ramp_agc; ///<
    RealT ramp_10;  ///<
    RealT ramp_30;  ///<
    RealT ramp_q;   ///<
    RealT apf;      ///<

    inline std::string str() const
    {
      std::stringstream ss;
      ss << std::setw(10) << bus 
         << std::setw(10) << Pg
         << std::setw(10) << Qg
         << std::setw(10) << Qmax
         << std::setw(10) << Qmin
         << std::setw(10) << Vg
         << std::setw(10) << mBase
         << std::setw(10) << status
         << std::setw(10) << Pmax
         << std::setw(10) << Pmin
         << std::setw(10) << Pc1
         << std::setw(10) << Pc2
         << std::setw(10) << Qc1min
         << std::setw(10) << Qc1max
         << std::setw(10) << Qc2min
         << std::setw(10) << Qc2max
         << std::setw(10) << ramp_agc
         << std::setw(10) << ramp_10
         << std::setw(10) << ramp_30
         << std::setw(10) << ramp_q
         << std::setw(10) << apf;
      ss << "\n";
      return ss.str();
    }
  };

  template <typename RealT = double, typename IdxT = int>
  struct BranchData
  {
    IdxT  fbus;   ///< "From" bus ID
    IdxT  tbus;   ///< "To" bus ID
    RealT r;      ///< Resistance (p.u.)
    RealT x;      ///< Reactance (p.u.)
    RealT b;      ///< Total line charging susceptance (p.u.)
    RealT rateA;  ///< MVA rating A (long term rating), 0=unlimited
    RealT rateB;  ///< MVA rating B (short term rating), 0=unlimited
    RealT rateC;  ///< MVA rating C (emergency rating), 0=unlimited
    RealT ratio;  ///< Transformer off nominal turns ratio
    RealT angle;  ///< Transformer phase shift angle [deg], positive ⇒ delay
    IdxT  status; ///< Initial service status: 1=in-service, 0=out-of-service
    RealT angmin; ///< Minimum anngle difference af - at [deg]
    RealT angmax; ///< Maximum anngle difference af - at [deg]

    inline std::string str() const
    {
      std::stringstream ss;
      ss << std::setw(10) << fbus
         << std::setw(10) << tbus
         << std::setw(10) << r
         << std::setw(10) << x
         << std::setw(10) << b
         << std::setw(10) << rateA
         << std::setw(10) << rateB
         << std::setw(10) << rateC
         << std::setw(10) << ratio
         << std::setw(10) << angle
         << std::setw(10) << status
         << std::setw(10) << angmin
         << std::setw(10) << angmax;
      ss << "\n";
      return ss.str();
    }
  };

  template <typename RealT = double, typename IdxT = int>
  struct GenCostData
  {
    IdxT kind;
    IdxT startup;
    IdxT shutdown;
    IdxT n;
    std::vector<RealT> rest;

    inline std::string str() const
    {
      std::stringstream ss;
      ss << std::setw(10) << kind
         << std::setw(10) << startup
         << std::setw(10) << shutdown
         << std::setw(10) << n;
      for (const auto& val : rest)
        ss << std::setw(10) << val;
      ss << "\n";
      return ss.str();
    }
  };

  template <typename RealT = double, typename IdxT = int>
  struct SystemModelData
  {
    using BusDataT = BusData<RealT, IdxT>;
    using GenDataT = GenData<RealT, IdxT>;
    using BranchDataT = BranchData<RealT, IdxT>;
    using GenCostDataT = GenCostData<RealT, IdxT>;
    using LoadDataT = LoadData<RealT, IdxT>;

    std::string version;
    IdxT baseMVA;
    std::vector<BusDataT> bus;
    std::vector<GenDataT> gen;
    std::vector<BranchDataT> branch;
    std::vector<GenCostDataT> gencost;
    std::vector<LoadDataT> load;

    // Not sure if these should be in this struct... Not all matpower files
    // I found contained them.
    //
    // std::string name;
    // std::vector<std::string> bus_name;

    inline std::string str() const
    {
      std::stringstream ss;
      ss << "Version: " << version << "\n"
         << "Base MVA: " << baseMVA << "\n";

      ss << "Bus:\n";
      for (const auto& v : bus)
        ss << bus.str() << "\n";

      ss << "Gen:\n";
      for (const auto& v : gen)
        ss << gen.str();

      ss << "Branch:\n";
      for (const auto& v : branch)
        ss << branch.str();

      ss << "GenCost:\n";
      for (const auto& v : gencost)
        ss << gencost.str();

      ss << "\n";

      return ss.str();
    }
  }; // struct SystemModelData

}  // namespace PowerSystemData
}  // namespace GridKit
