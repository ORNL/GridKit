#include <array>
#include <cmath>
#include <span>
#include <tuple>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKitTest
{
  using scalar_type = double;
  using index_type  = size_t;

  bool close(scalar_type value, scalar_type ref, scalar_type tol = 1.0e-10)
  {
    return std::abs(value - ref) <= tol * (1.0 + std::abs(ref));
  }

  class PlumbingElement : public GridKit::EMT::Element<scalar_type, index_type>
  {
  public:
    using BaseT          = GridKit::EMT::Element<scalar_type, index_type>;
    using SignalT        = typename BaseT::SignalT;
    using SignalStorageT = typename BaseT::SignalStorageT;

    index_type size() const override
    {
      return 2;
    }

    std::span<const SignalT> inputs() const override
    {
      return std::span<const SignalT>(inputs_.data(), inputs_.size());
    }

    int initialize() override
    {
      return 0;
    }

    void tagDifferentiable(std::vector<bool>& tag) const override
    {
      for (index_type i = base_; i < end(); ++i)
      {
        tag[i] = true;
      }
    }

    template <typename V>
    void evaluateResidualKernel(V* y, V* yp, V* u, V* f) const
    {
      auto state_input      = input(u, inputs_[0]);
      auto derivative_input = input(u, inputs_[1]);

      f[0] = -yp[0] + y[0] * y[1] + state_input[0] + derivative_input[0];
      f[1] = -y[1] + 2.0 * y[0] + 3.0 * state_input[0] - 4.0 * derivative_input[0];
    }

    int evaluateResidual() override
    {
      evaluateResidualKernel<scalar_type>(y_, yp_, nullptr, f_);
      return 0;
    }

    int evaluateJacobian() override
    {
      return this->template evaluateInputDifferentialJacobian<PlumbingElement>();
    }

  private:
    std::array<SignalT, 2> inputs_{
        SignalT{0, 1, 1, SignalStorageT::State},
        SignalT{1, 1, 1, SignalStorageT::Derivative}};
  };

  scalar_type cooValue(PlumbingElement::MatrixT& J, index_type row, index_type col)
  {
    auto [rows, cols, vals] = J.getEntryCopies();
    scalar_type value       = 0.0;
    for (size_t i = 0; i < vals.size(); ++i)
    {
      if (rows[i] == row && cols[i] == col)
      {
        value += vals[i];
      }
    }
    return value;
  }

  GridKit::Testing::TestOutcome sparse_jet_owned_and_borrowed_inputs()
  {
    using GridKit::Testing::TestStatus;

    PlumbingElement element;

    std::array<scalar_type, 4> y{11.0, 0.0, 4.0, 5.0};
    std::array<scalar_type, 4> yp{0.0, 13.0, 1.5, 0.0};
    std::array<scalar_type, 4> f{};
    std::array<index_type, 4>  index{0, 1, 2, 3};

    element.bind(y.data(), yp.data(), f.data(), index.data(), 2);
    element.setCoordinate(60.0, 7.0);

    TestStatus success = true;

    element.evaluateResidual();
    success *= close(f[2], -1.5 + 4.0 * 5.0 + 11.0 + 13.0);
    success *= close(f[3], -5.0 + 2.0 * 4.0 + 3.0 * 11.0 - 4.0 * 13.0);

    element.evaluateJacobian();
    auto& J = element.jacobian();

    success *= (J.nnz() > 0);
    success *= close(cooValue(J, 2, 2), -2.0);
    success *= close(cooValue(J, 2, 3), 4.0);
    success *= close(cooValue(J, 2, 0), 1.0);
    success *= close(cooValue(J, 2, 1), 7.0);
    success *= close(cooValue(J, 3, 2), 2.0);
    success *= close(cooValue(J, 3, 3), -1.0);
    success *= close(cooValue(J, 3, 0), 3.0);
    success *= close(cooValue(J, 3, 1), -28.0);

    return success.report(__func__);
  }
} // namespace GridKitTest

int main()
{
  GridKit::Testing::TestingResults result;
  result += GridKitTest::sparse_jet_owned_and_borrowed_inputs();
  return result.summary();
}
