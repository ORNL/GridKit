/**
 * @file runModelTests.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Runs the EMT model tests.
 */

#include "ModelTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::ModelTests     test;

  result += test.loadZ();
  result += test.loadRLC();
  result += test.invalidLoadRLC();
  result += test.conjugatePoles();
  result += test.invalidVectorFit();
  result += test.invalidInput();
  result += test.absoluteTolerances();
  result += test.monitors();
  result += test.operatorDevice();
  result += test.heldInput();
  result += test.variableReference();
  result += test.outputReference();
  result += test.pmu();
  result += test.pmuReference();
  result += test.layout();
  result += test.voltageSource();
  result += test.lineLumped();
  result += test.bundledLineLumped();
  result += test.invalidLineLumped();
  result += test.closedSwitch();
  result += test.openSwitch();
  result += test.switchWaitsForZero();
  result += test.sharedCommand();
  result += test.machine();
  result += test.machineIronLoss();
  result += test.transformer();
  result += test.machineCurrentReference();
  result += test.pll();
  result += test.pllFrequencyReference();
  result += test.regca();
  result += test.park();
  result += test.tgov1();
  result += test.gastPti();
  result += test.hygov();
  result += test.sexsPti();
  result += test.ieeet1();
  result += test.esdc1a();
  result += test.ieeest();
  result += test.reecb();
  result += test.repca();
  result += test.machineControls();

  return result.summary();
}
