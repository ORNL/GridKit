

#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <map>

/**
 * @brief A very basic hypergraph setup for circuit representation.
 * This forms the hypergraph as a bipartite graph. Doesn't allow
 * removing. Can only grab sets of connections to nodes
 *
 * @todo should replace with something better and more efficent.
 * Should replace with a libraries setup instead. This would allow
 * fast and easy partitioning of circuits
 * 
 * @todo This is to replace inserting vector size for allocating PowerElectronicsModel
 *
 * @todo should replace N and E with Node and Component classes respectively.
 *
 * @note Tested but currently not used in the rest of the code.
 *
 * @tparam IdxT
 * @tparam Label
 */
template <typename N, typename E>
class CircuitGraph
{
private:
	std::set<N> hypernodes;
	std::set<E> hyperedges;
	std::map<E, std::set<N>> edgestonodes;

public:
	CircuitGraph();
	~CircuitGraph();
	bool addHyperEdge(E he);
	bool addHyperNode(N hn);
	bool addConnection(N hn, E he);
	std::set<N> getHyperEdgeConnections(E he);
	size_t amountHyperNodes();
	size_t amountHyperEdges();
	void printBiPartiteGraph(bool verbose = false);
};

template <typename N, typename E>
CircuitGraph<N, E>::CircuitGraph()
{
}

template <typename N, typename E>
CircuitGraph<N, E>::~CircuitGraph()
{
}

template <typename N, typename E>
bool CircuitGraph<N, E>::addHyperNode(N hn)
{
	return this->hypernodes.insert(hn).second;
}

template <typename N, typename E>
bool CircuitGraph<N, E>::addHyperEdge(E he)
{
	return this->hyperedges.insert(he).second;
}

template <typename N, typename E>
bool CircuitGraph<N, E>::addConnection(N hn, E he)
{
	if (this->hyperedges.count(he) == 0 || this->hypernodes.count(hn) == 0)
	{
		return false;
	}
	return this->edgestonodes[he].insert(hn).second;
}

template <typename N, typename E>
std::set<N> CircuitGraph<N, E>::getHyperEdgeConnections(E he)
{
	return this->edgestonodes[he];
}

template <typename N, typename E>
size_t CircuitGraph<N, E>::amountHyperNodes()
{
	return this->hypernodes.size();
}

template <typename N, typename E>
size_t CircuitGraph<N, E>::amountHyperEdges()
{
	return this->hyperedges.size();
}

/**
 * @brief Printing
 *
 * @todo need to add verbose printing for connections display
 *
 * @tparam IdxT
 * @param verbose
 */

template <typename N, typename E>
void CircuitGraph<N, E>::printBiPartiteGraph(bool verbose)
{

	std::cout << "Amount of HyperNodes: " << this->amountHyperNodes() << std::endl;
	std::cout << "Amount of HyperEdges: " << this->amountHyperEdges() << std::endl;
	std::cout << "Connections per Edge:" << std::endl;
	for (auto i : this->edgestonodes)
	{
		std::cout << i.first << " : {";
		for (auto j : i.second)
		{
			std::cout << j << ", ";
		}
		std::cout << "}\n";
	}
}
