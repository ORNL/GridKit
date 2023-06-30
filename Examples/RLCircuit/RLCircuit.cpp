

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
#include <Node.hpp>


int main(int argc, char const *argv[])
{
	//Basic Circuit Setup
	//@todo will want to eventually put components in the hyper graph instead of indicies

	 std::unique_ptr<CircuitGraph<size_t, size_t>> cirg(new CircuitGraph<size_t,size_t>());
	for (size_t i = 0; i < 3; i++)
	{
		cirg->addHyperEdge(i);
	}
	for (size_t i = 0; i < 5; i++)
	{
		cirg->addHyperNode(i);
	}

	//Create Connections of Nodes to edges sets
	//External nodes
	cirg->addConnection(0, 0);
	cirg->addConnection(0, 2);
	cirg->addConnection(2, 0);
	cirg->addConnection(2, 1);
	cirg->addConnection(3, 1);
	cirg->addConnection(3, 2);

	//Internal nodes
	cirg->addConnection(1,0);
	cirg->addConnection(4,2);


	cirg->printBiPartiteGraph();
	

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
	sysmodel->allocate(*cirg);

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


	return 0;
}
