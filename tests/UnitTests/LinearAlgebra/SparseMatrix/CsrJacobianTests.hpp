#include <sstream>

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Utilities/TestHelpers.hpp>
#include <GridKit/Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class CSRJacobianTests
    {
    public:
      /**
       * @brief Test the COO jacobian vs CSR jacobian of a component, and ensure that they agree.
       * @param comp_name The name of the component being checked, for reporting purposes.
       * @param args The arguments used to construct the component
       */
      template <template <class S, typename I> typename Component, class... Args>
      TestOutcome testCooVsCsrJacobian(std::string comp_name, Args&&... args)
      {
        using namespace Model;
        using namespace LinearAlgebra;

        using COO_Matrix  = COO_Matrix<ScalarT, IdxT>;
        using CSRJacobian = Evaluator<ScalarT, IdxT>::CSRJacobian;

        TestStatus success = true;

        Component<ScalarT, IdxT> comp(args...);

        if (!comp.hasJacobian())
        {
          std::cerr << "Testing Jacobian of " << comp_name << ", but none was available." << std::endl;
          success = false;
        }

        if (!comp.hasCSRJacobian())
        {
          std::cerr << "Testing CSR Jacobian of " << comp_name << ", but none was available." << std::endl;
          success = false;
        }

        comp.evaluateJacobian();
        comp.evaluateCSRJacobian();

        COO_Matrix&  coo_jac = comp.getJacobian();
        CSRJacobian& csr_jac = comp.getCSRJacobian();

        success *= compare(coo_jac, csr_jac);

        std::stringstream test_name;
        test_name << __func__ << " <" << comp_name << ">";
        return success.report(test_name.str().c_str());
      }

    private:
      bool compare(LinearAlgebra::COO_Matrix<ScalarT, IdxT>& original_mat, LinearAlgebra::CSRMatrix<ScalarT, IdxT>& new_mat)
      {
        bool comparison = true;

        auto [rows, cols, vals] = original_mat.getEntries();

        comparison = comparison
                     && std::get<0>(original_mat.getDimensions()) == new_mat.numRows()
                     && std::get<1>(original_mat.getDimensions()) == new_mat.numCols()
                     && vals.size() == new_mat.numNonZero();

        if (!comparison)
          return comparison;

        ScalarT target;
        // Test all elements - make sure they match
        for (size_t i = 0; i < new_mat.numRows(); i++)
        {
          for (IdxT j = 0; j < new_mat.numCols(); j++)
          {
            target = 0;
            for (size_t k = 0; k < vals.size(); k++)
            {
              if (rows[k] == i && cols[k] == j)
              {
                target = vals[k];
                break;
              }
            }

            comparison = comparison && new_mat.valueAt(i, j) == target;
          }
        }

        return comparison;
      }
    };
  } // namespace Testing
} // namespace GridKit
