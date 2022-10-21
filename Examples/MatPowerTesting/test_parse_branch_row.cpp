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

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.00281	0.0281	0.00712	400	400	400	0	0	1	-360	360;
	1	4	0.00304	0.0304	0.00658	0	0	0	0	0	1	-360	360;
	1	5	0.00064	0.0064	0.03126	0	0	0	0	0	1	-360	360;
	2	3	0.00108	0.0108	0.01852	0	0	0	0	0	1	-360	360;
	3	4	0.00297	0.0297	0.00674	0	0	0	0	0	1	-360	360;
	4	5	0.00297	0.0297	0.00674	240	240	240	0	0	1	-360	360;
];


)"};

} // namespace

int main(int argc, char **argv) {
  int fail = 0;
  std::vector<BranchRow<RealT, IdxT>> branch_answer{
      {1, 2, 0.00281, 0.0281, 0.00712, 400, 400, 400, 0, 0, 1, -360, 360},
      {1, 4, 0.00304, 0.0304, 0.00658, 0, 0, 0, 0, 0, 1, -360, 360},
      {1, 5, 0.00064, 0.0064, 0.03126, 0, 0, 0, 0, 0, 1, -360, 360},
      {2, 3, 0.00108, 0.0108, 0.01852, 0, 0, 0, 0, 0, 1, -360, 360},
      {3, 4, 0.00297, 0.0297, 0.00674, 0, 0, 0, 0, 0, 1, -360, 360},
      {4, 5, 0.00297, 0.0297, 0.00674, 240, 240, 240, 0, 0, 1, -360, 360},
  };
  MatPower<RealT, IdxT> mp;

  {
    std::istringstream iss(matpower_data);
    GridKit::readMatPower(mp, iss);
    if (!isEqual(mp.branch, branch_answer))
      fail++;
    std::cout << "After reading the branch component, fail == " << fail << "\n";
  }

  std::cout << "Tests " << (fail ? "FAILED" : "PASSED") << "\n";

  return fail;
}
