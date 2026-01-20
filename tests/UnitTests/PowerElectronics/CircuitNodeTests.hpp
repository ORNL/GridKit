
#include <iomanip>
#include <iostream>

#include <GridKit/Model/PowerElectronics/CircuitNode.hpp>
#include <GridKit/Utilities/TestHelpers.hpp>
#include <GridKit/Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class NodeTests
    {
    public:
      NodeTests()  = default;
      ~NodeTests() = default;

      /// Constructor, allocation, and initialization checks
      TestOutcome constructor()
      {
        TestStatus success = true;

        ScalarT V{1.0};

        CircuitNode<ScalarT, IdxT>* node = nullptr;

        // Default construct
        node = new CircuitNode<ScalarT, IdxT>();
        node->allocate();
        node->initialize();
        success *= isEqual(node->V(), static_cast<ScalarT>(0));
        success *= isEqual(node->I(), static_cast<ScalarT>(0));
        delete node;

        // Construct with initial voltage
        node = new CircuitNode<ScalarT, IdxT>(V);
        node->allocate();
        node->initialize();
        success *= isEqual(node->V(), V);
        success *= isEqual(node->I(), static_cast<ScalarT>(0));
        delete node;

        node = nullptr;

        return success.report(__func__);
      }

      /// Accessor method tests
      TestOutcome residual()
      {
        TestStatus success = true;

        ScalarT V{1.0};
        ScalarT I{1.0};

        CircuitNode<ScalarT, IdxT> node(V);
        node.allocate();
        node.initialize();
        success *= isEqual(node.V(), V);

        node.I()  = I;
        success  *= isEqual(node.I(), I);

        node.evaluateResidual();
        success *= isEqual(node.I(), 0.0);

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
