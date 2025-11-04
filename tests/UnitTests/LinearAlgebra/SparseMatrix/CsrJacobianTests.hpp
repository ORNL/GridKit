#include <format>
#include <sstream>

#include "LinearAlgebra/SparseMatrix/CsrMatrix.hpp"
#include "Model/Evaluator.hpp"
#include "Utilities/TestHelpers.hpp"
#include "Utilities/Testing.hpp"

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class CsrJacobianTests
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
        using CsrJacobian = Evaluator<ScalarT, IdxT>::CsrJacobian;

        TestStatus success = true;

        Component<ScalarT, IdxT> comp(args...);

        if (!comp.hasJacobian())
        {
          std::cerr << "Testing Jacobian of " << comp_name << ", but none was available." << std::endl;
          success = false;
        }

        if (!comp.hasCsrJacobian())
        {
          std::cerr << "Testing CSR Jacobian of " << comp_name << ", but none was available." << std::endl;
          success = false;
        }

        comp.evaluateJacobian();
        comp.evaluateCsrJacobian();

        COO_Matrix&  coo_jac = comp.getJacobian();
        CsrJacobian& csr_jac = comp.getCsrJacobian();

        success *= compare(coo_jac, csr_jac);

        std::stringstream test_name;
        test_name << __func__ << " <" << comp_name << ">";
        return success.report(test_name.str().c_str());
      }

      TestOutcome testSystemCooVsCsrJacobian(std::string model_name, Model::Evaluator<ScalarT, IdxT>& model)
      {
        using namespace Model;
        using namespace LinearAlgebra;

        using COO_Matrix  = COO_Matrix<ScalarT, IdxT>;
        using CsrJacobian = Evaluator<ScalarT, IdxT>::CsrJacobian;

        TestStatus success = true;

        model.evaluateJacobian();
        model.evaluateCsrJacobian();

        COO_Matrix&  coo_jac = model.getJacobian();
        CsrJacobian& csr_jac = model.getCsrJacobian();

        success *= compare(coo_jac, csr_jac);

        std::stringstream test_name;
        test_name << __func__ << " <" << model_name << ">";
        return success.report(test_name.str().c_str());
      }

    private:
      bool compare(LinearAlgebra::COO_Matrix<ScalarT, IdxT>& original_mat, LinearAlgebra::CsrMatrix<ScalarT, IdxT>& new_mat)
      {
        bool comparison = true;

        auto [rows, cols, vals] = original_mat.getEntries();

        comparison = comparison
                     && std::get<0>(original_mat.getDimensions()) == new_mat.numRows()
                     && std::get<1>(original_mat.getDimensions()) == new_mat.numCols();

        if (!comparison)
        {
          std::cerr << std::format("Matrix dimensions do not align:\n{},{} v.s.\n{},{}\n",
                                   std::get<0>(original_mat.getDimensions()),
                                   std::get<1>(original_mat.getDimensions()),
                                   new_mat.numRows(),
                                   new_mat.numCols());

          original_mat.printMatrix();
          std::cerr << new_mat.printNonzeroElements() << '\n';
          return comparison;
        }

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
