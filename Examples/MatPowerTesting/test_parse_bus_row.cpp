#include <FileIO.hpp>
#include <PowerSystemData.hpp>
#include <Testing.hpp>
#include <iostream>

using namespace GridKit;
using namespace GridKit::Testing;
using namespace GridKit::PowerSystemData;

namespace {

using IdxT = int;
using RealT = double;

static const std::string matpower_data{
    R"(
%% bus data
mpc.bus = [
	1	2	0	0	0	0	1	1	0	230	1	1.1	0.0;
	2	1	300	98.61	0	0	1	1	0	230	1	1.1	0.0;
	3	2	300	98.61	0	0	1	1	0	230	1	1.1	0.0;
	4	3	400	131.47	0	0	1	1	0	230	1	1.1	0.0;
	5	2	0	0	0	0	1	1	0	230	1	1.1	0.9;
];
)"};

}  // namespace

int main(int argc, char** argv) {
  int fail = 0;
  std::vector<BusData<RealT, IdxT>> bus_answer{
      {1, 2, 0, 0, 1, 1, 0, 230, 1, 1.1, 0.0},
      {2, 1, 0, 0, 1, 1, 0, 230, 1, 1.1, 0.0},
      {3, 2, 0, 0, 1, 1, 0, 230, 1, 1.1, 0.0},
      {4, 3, 0, 0, 1, 1, 0, 230, 1, 1.1, 0.0},
      {5, 2, 0, 0, 1, 1, 0, 230, 1, 1.1, 0.9},
  };

  std::vector<LoadData<RealT, IdxT>> load_answer{
      {1,   0,      0},
      {2, 300,  98.61},
      {3, 300,  98.61},
      {4, 400, 131.47},
      {5,   0,      0},
  };

  {
    SystemModelData<RealT, IdxT> mp;
    std::istringstream iss(matpower_data);
    GridKit::readMatPower(mp, iss);
    if (!isEqual(mp.bus, bus_answer)) fail++;
    std::cout << "After reading the bus component, fail == " << fail << "\n";
  }

  {
    SystemModelData<RealT, IdxT> mp;
    std::istringstream iss(matpower_data);
    GridKit::readMatPower(mp, iss);
    if (!isEqual(mp.load, load_answer)) fail++;
    std::cout << "After reading the bus component, fail == " << fail << "\n";
  }

  std::cout << "Tests " << (fail?"FAILED":"PASSED") << "\n";

  return fail;
}
