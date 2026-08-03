#include <array>
#include <cmath>

#include <GridKit/LinearAlgebra/DenseExpression.hpp>
#include <GridKit/LinearAlgebra/SparseJet.hpp>
#include <GridKit/Model/EMT/Signal.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using scalar_type = double;
  using index_type  = size_t;
  using SliceT      = GridKit::EMT::Slice<scalar_type, index_type>;
  using JetT        = GridKit::LinearAlgebra::SparseJet<scalar_type, index_type>;
  using JetSliceT   = GridKit::EMT::Slice<JetT, index_type>;

  bool close(scalar_type value, scalar_type ref, scalar_type tol = 1.0e-12)
  {
    return std::abs(value - ref) <= tol * (1.0 + std::abs(ref));
  }

  GridKit::Testing::TestOutcome dense_expression_smoke()
  {
    using GridKit::Testing::TestStatus;
    using namespace GridKit::LinearAlgebra;

    std::array<index_type, 4>  matrix_index{0, 1, 2, 3};
    std::array<index_type, 2>  vector_index{0, 1};
    std::array<scalar_type, 4> a_values{1.0, 2.0, 3.0, 4.0};
    std::array<scalar_type, 4> b_values{5.0, 6.0, 7.0, 8.0};
    std::array<scalar_type, 2> v_values{2.0, 3.0};
    std::array<scalar_type, 4> out_values{};
    std::array<scalar_type, 2> vec_out_values{};

    SliceT A{a_values.data(), matrix_index.data(), 2, 2};
    SliceT B{b_values.data(), matrix_index.data(), 2, 2};
    SliceT V{v_values.data(), vector_index.data(), 2, 1};
    SliceT Out{out_values.data(), matrix_index.data(), 2, 2};
    SliceT VecOut{vec_out_values.data(), vector_index.data(), 2, 1};

    TestStatus success = true;

    equation(Out)  = -A + B;
    success       *= close(Out(0, 0), 4.0);
    success       *= close(Out(1, 1), 4.0);

    equation(Out)  = A * B + 2.0 * identity<scalar_type, index_type>(2);
    success       *= close(Out(0, 0), 21.0);
    success       *= close(Out(0, 1), 22.0);
    success       *= close(Out(1, 0), 43.0);
    success       *= close(Out(1, 1), 52.0);

    equation(Out)  = diag(V) * A - hadamard(A, B) / 2.0;
    success       *= close(Out(0, 0), -0.5);
    success       *= close(Out(0, 1), -2.0);
    success       *= close(Out(1, 0), -1.5);
    success       *= close(Out(1, 1), -4.0);

    auto M            = matrix<scalar_type, index_type>(2,
                                             2,
                                             [](index_type i, index_type j)
                                             {
                                               return (i == j) ? 4.0 : 1.0;
                                             });
    auto v            = vector<scalar_type, index_type>(2,
                                             [](index_type i)
                                             {
                                               return static_cast<scalar_type>(i + 1);
                                             });
    equation(VecOut)  = M * v;
    success          *= close(VecOut[0], 6.0);
    success          *= close(VecOut[1], 9.0);

    zero(Out);
    success *= close(Out(0, 0), 0.0);
    success *= close(Out(1, 1), 0.0);

    return success.report(__func__);
  }

  scalar_type tangent_value(const JetT& jet, index_type col)
  {
    scalar_type value = 0.0;
    for (const auto& [entry_col, entry_value] : jet.tangent)
    {
      if (entry_col == col)
      {
        value += entry_value;
      }
    }
    return value;
  }

  GridKit::Testing::TestOutcome sparse_jet_smoke()
  {
    using GridKit::Testing::TestStatus;
    using namespace GridKit::LinearAlgebra;

    std::array<index_type, 4> matrix_index{0, 1, 2, 3};
    std::array<JetT, 4>       p_values{
        JetT::variable(1.0, 10),
        JetT::variable(2.0, 11),
        JetT::variable(3.0, 12),
        JetT::variable(4.0, 13)};
    std::array<JetT, 4> c_values{
        JetT::variable(5.0, 20),
        JetT::variable(6.0, 21),
        JetT::variable(7.0, 22),
        JetT::variable(8.0, 23)};
    std::array<JetT, 4> out_values{};

    JetSliceT P{p_values.data(), matrix_index.data(), 2, 2};
    JetSliceT C{c_values.data(), matrix_index.data(), 2, 2};
    JetSliceT Out{out_values.data(), matrix_index.data(), 2, 2};

    TestStatus success = true;

    equation(Out) = P * C - identity(P);

    success *= close(Out(0, 0).value, 18.0);
    success *= close(tangent_value(Out(0, 0), 10), 5.0);
    success *= close(tangent_value(Out(0, 0), 11), 7.0);
    success *= close(tangent_value(Out(0, 0), 20), 1.0);
    success *= close(tangent_value(Out(0, 0), 22), 2.0);
    success *= close(tangent_value(Out(0, 0), 13), 0.0);

    auto x   = JetT::variable(2.0, 31);
    auto y   = x + x;
    success *= close(y.value, 4.0);
    success *= close(tangent_value(y, 31), 2.0);

    auto masked  = x * 0.0;
    success     *= close(masked.value, 0.0);
    success     *= (masked.tangent.size() == 1);
    success     *= close(tangent_value(masked, 31), 0.0);

    using GridKit::LinearAlgebra::atan2;
    using GridKit::LinearAlgebra::log;
    auto z   = log(x * x) + atan2(x, JetT(3.0));
    success *= close(z.value, std::log(4.0) + std::atan2(2.0, 3.0));
    success *= close(tangent_value(z, 31), 1.0 + 3.0 / 13.0);

    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += dense_expression_smoke();
  result += sparse_jet_smoke();
  return result.summary();
}
