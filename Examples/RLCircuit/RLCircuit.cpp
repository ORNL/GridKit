

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>

#include <ComponentLib/PowerElectronicsComponents/Capacitor/Capacitor.hpp>
#include <ComponentLib/PowerElectronicsComponents/Inductor/Inductor.hpp>
#include <ComponentLib/PowerElectronicsComponents/Resistor/Resistor.hpp>
#include <ComponentLib/PowerElectronicsComponents/VoltageSource/VoltageSource.hpp>

#include <PowerElectronicsModel.hpp>


int main(int argc, char const *argv[])
{
	

	//Create circuit model
	ModelLib::PowerElectronicsModel<double, size_t>* sysmodel = new ModelLib::PowerElectronicsModel<double, size_t>();

	size_t idoff = 0;

	//inductor
	ModelLib::Inductor<double, size_t>* induct = new ModelLib::Inductor(idoff,0.1);
	//Form index to node uid realations
	induct->setExternalConnectionNodes(0,2);
	induct->setExternalConnectionNodes(2,1);
	induct->setExternalConnectionNodes(1,0);
	sysmodel->addComponent(induct);


	//resistor
	idoff++;
	ModelLib::Resistor<double, size_t>* resis = new ModelLib::Resistor(idoff, 1.0);
	//Form index to node uid realations
	resis->setExternalConnectionNodes(0,3);
	resis->setExternalConnectionNodes(1,2);
	sysmodel->addComponent(resis);

	//voltage source
	idoff++;
	ModelLib::VoltageSource<double, size_t>* vsource = new ModelLib::VoltageSource(idoff, 0.1);
	//Form index to node uid realations
	vsource->setExternalConnectionNodes(0,0);
	vsource->setExternalConnectionNodes(2,4);
	vsource->setExternalConnectionNodes(1,3);
	sysmodel->addComponent(vsource);


	//Allocate with graph
	sysmodel->allocate(5);

	std::cout << sysmodel->y().size() << std::endl;

	//Create Intial points
	sysmodel->y()[0] = 1.0;
	sysmodel->y()[1] = 1.0;
	sysmodel->y()[2] = 1.0;
	sysmodel->y()[3] = 1.0;
	sysmodel->y()[4] = 1.0;

	sysmodel->yp()[0] = 1.0;
	sysmodel->yp()[1] = 1.0;
	sysmodel->yp()[2] = 1.0;
	sysmodel->yp()[3] = 1.0;
	sysmodel->yp()[4] = 1.0;

	
	sysmodel->initialize();

	sysmodel->evaluateResidual();

	std::cout << "Output: {";
	for (double i : sysmodel->getResidual())
	{
		std::cout << i << ", ";
	}
	std::cout << "}\n";

	induct->updateTime(0.0, 1.0);
	induct->evaluateJacobian();
	induct->getJacobian().printMatrix(true);


	return 0;
}
