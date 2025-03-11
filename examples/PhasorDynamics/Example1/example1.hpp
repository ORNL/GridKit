#include<stdio.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<time.h>

//#include <sundials_core.h>
#include <idas/idas.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunlinsol/sunlinsol_dense.h>

#include "Model/PhasorDynamics/SystemModel.hpp"
#include "Model/PhasorDynamics/Bus/Bus.hpp"
#include "Model/PhasorDynamics/Bus/BusInfinite.hpp"
#include "Model/PhasorDynamics/Branch/Branch.hpp"
#include "Model/PhasorDynamics/BusFault/BusFault.hpp"
#include "Model/PhasorDynamics/SynchronousMachine/GENROUwS/GENROU.hpp"


#include "Solver/Dynamic/Ida.hpp"

#include "Model/PhasorDynamics/Bus/Bus.cpp"
#include "Model/PhasorDynamics/Bus/BusInfinite.cpp"
#include "Model/PhasorDynamics/Branch/Branch.cpp"
#include "Solver/Dynamic/Ida.cpp"



int main();
