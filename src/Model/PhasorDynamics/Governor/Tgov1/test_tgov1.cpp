#include <iostream>

#include "Tgov1.hpp"
#include "Tgov1Data.hpp"

int main()
{
  using ScalarT = double;
  using IdxT    = size_t;

  // Test machine pointer for testing, can be nullptr or mock object
  GridKit::PhasorDynamics::Genrou<ScalarT, IdxT>* dummy_machine = nullptr;

  GridKit::PhasorDynamics::Governor::Tgov1Data<ScalarT, IdxT> data;

  // Set parameters by inserting values into the parameters map with correct enum keys
  data.parameters[GridKit::PhasorDynamics::Governor::Tgov1Parameters::R]     = ScalarT(0.04);
  data.parameters[GridKit::PhasorDynamics::Governor::Tgov1Parameters::Pvmin] = ScalarT(0.0);
  data.parameters[GridKit::PhasorDynamics::Governor::Tgov1Parameters::Pvmax] = ScalarT(1.2);
  data.parameters[GridKit::PhasorDynamics::Governor::Tgov1Parameters::T1]    = ScalarT(0.5);
  data.parameters[GridKit::PhasorDynamics::Governor::Tgov1Parameters::T2]    = ScalarT(3.0);
  data.parameters[GridKit::PhasorDynamics::Governor::Tgov1Parameters::T3]    = ScalarT(6.0);
  data.parameters[GridKit::PhasorDynamics::Governor::Tgov1Parameters::Dt]    = ScalarT(0.01);

  // Construct Tgov1 with this data
  GridKit::PhasorDynamics::Governor::Tgov1<ScalarT, IdxT> tgov1(dummy_machine, data);

  std::cout << "Tgov1 constructed successfully." << std::endl;

  return 0;
}
