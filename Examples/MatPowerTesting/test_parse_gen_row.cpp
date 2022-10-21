#include <FileIO.hpp>
#include <MatPowerUtils.hpp>
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
%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	40	0	30	-30	1	100	1	40	0	0	0	0	0	0	0	0	0	0	0	0;
	1	170	0	127.5	-127.5	1	100	1	170	0	0	0	0	0	0	0	0	0	0	0	0;
	3	323.49	0	390	-390	1	100	1	520	0	0	0	0	0	0	0	0	0	0	0	0;
	4	0	0	150	-150	1	100	1	200	0	0	0	0	0	0	0	0	0	0	0	0;
	5	466.51	0	450	-450	1	100	1	600	0	0	0	0	0	0	0	0	0	0	0	0;
];

)"};

} // namespace

int main(int argc, char **argv) {
  int fail = 0;
  std::vector<GenRow<RealT, IdxT>> gen_answer{
      {1, 40, 0, 30, -30, 1, 100, 1, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1, 170, 0, 127.5, -127.5, 1, 100, 1, 170, 0, 0,
       0, 0,   0, 0,     0,      0, 0,   0, 0,   0},
      {3, 323.49, 0, 390, -390, 1, 100, 1, 520, 0, 0,
       0, 0,      0, 0,   0,    0, 0,   0, 0,   0},
      {4, 0, 0, 150, -150, 1, 100, 1, 200, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {5, 466.51, 0, 450, -450, 1, 100, 1, 600, 0, 0,
       0, 0,      0, 0,   0,    0, 0,   0, 0,   0}};
  MatPower<RealT, IdxT> mp;

  {
    std::istringstream iss(matpower_data);
    GridKit::readMatPower(mp, iss);
    if (!isEqual(mp.gen, gen_answer))
      fail++;
    std::cout << "After reading the gen component, fail == " << fail << "\n";
  }

  std::cout << "Tests " << (fail ? "FAILED" : "PASSED") << "\n";

  return fail;
}
