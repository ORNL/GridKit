#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>

#include <IpIpoptApplication.hpp>

#include <GridKit/Solver/Optimization/Ipopt/IpoptAdapter.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "AffineModel.hpp"
#include <Model.hpp>

namespace GridKit
{
  namespace Testing
  {
    namespace Optimization
    {
      template <class scalar_type, typename index_type>
      class MissingJacobianModel : public Model<scalar_type, index_type>
      {
      public:
        bool hasJacobian() override
        {
          return false;
        }
      };

      template <class scalar_type, typename index_type>
      class MissingHessianModel : public Model<scalar_type, index_type>
      {
      public:
        bool hasHessian() override
        {
          return false;
        }
      };

      template <class scalar_type, typename index_type>
      class IpoptTests
      {
      public:
        using ScalarT  = scalar_type;
        using IdxT     = index_type;
        using ModelT   = Model<ScalarT, IdxT>;
        using RealT    = typename ModelT::RealT;
        using AdapterT = GridKit::Optimization::IpoptAdapter<ScalarT, IdxT>;

        TestOutcome exactIpoptCallbacks()
        {
          TestStatus success = true;
          ModelT     model;
          success *= model.allocate() == 0;
          success *= model.initialize() == 0;

          AdapterT adapter(model);

          Ipopt::Index                variable_count{0};
          Ipopt::Index                constraint_count{0};
          Ipopt::Index                jacobian_nonzeros{0};
          Ipopt::Index                hessian_nonzeros{0};
          Ipopt::TNLP::IndexStyleEnum index_style{};
          success *= adapter.get_nlp_info(variable_count,
                                          constraint_count,
                                          jacobian_nonzeros,
                                          hessian_nonzeros,
                                          index_style);
          success *= variable_count == 3;
          success *= constraint_count == 2;
          success *= jacobian_nonzeros == 5;
          success *= hessian_nonzeros == 5;
          success *= index_style == Ipopt::TNLP::C_STYLE;

          std::array<Ipopt::Number, 3> variable_lower_bounds{};
          std::array<Ipopt::Number, 3> variable_upper_bounds{};
          std::array<Ipopt::Number, 2> constraint_lower_bounds{};
          std::array<Ipopt::Number, 2> constraint_upper_bounds{};
          success *= adapter.get_bounds_info(variable_count,
                                             variable_lower_bounds.data(),
                                             variable_upper_bounds.data(),
                                             constraint_count,
                                             constraint_lower_bounds.data(),
                                             constraint_upper_bounds.data());
          success *= checkEntries(variable_lower_bounds,
                                  std::array<Ipopt::Number, 3>{{0.6, 0.6, 0.6}});
          success *= checkEntries(variable_upper_bounds,
                                  std::array<Ipopt::Number, 3>{{1.4, 1.4, 1.4}});
          success *= checkEntries(constraint_lower_bounds,
                                  std::array<Ipopt::Number, 2>{{0.0, 0.0}});
          success *= checkEntries(constraint_upper_bounds,
                                  std::array<Ipopt::Number, 2>{{0.0, 0.0}});

          std::array<Ipopt::Number, 3> starting_point{};
          success *= adapter.get_starting_point(variable_count,
                                                true,
                                                starting_point.data(),
                                                false,
                                                nullptr,
                                                nullptr,
                                                constraint_count,
                                                false,
                                                nullptr);
          success *= checkEntries(starting_point,
                                  std::array<Ipopt::Number, 3>{{0.8, 1.2, 0.9}});

          std::array<Ipopt::Index, 5> jacobian_rows{};
          std::array<Ipopt::Index, 5> jacobian_columns{};
          success *= adapter.eval_jac_g(variable_count,
                                        nullptr,
                                        false,
                                        constraint_count,
                                        jacobian_nonzeros,
                                        jacobian_rows.data(),
                                        jacobian_columns.data(),
                                        nullptr);
          success *= checkEntries(jacobian_rows,
                                  std::array<Ipopt::Index, 5>{{0, 0, 0, 1, 1}});
          success *= checkEntries(jacobian_columns,
                                  std::array<Ipopt::Index, 5>{{0, 1, 2, 1, 2}});

          std::array<Ipopt::Index, 5> hessian_rows{};
          std::array<Ipopt::Index, 5> hessian_columns{};
          success *= adapter.eval_h(variable_count,
                                    nullptr,
                                    false,
                                    1.0,
                                    constraint_count,
                                    nullptr,
                                    false,
                                    hessian_nonzeros,
                                    hessian_rows.data(),
                                    hessian_columns.data(),
                                    nullptr);
          success *= checkEntries(hessian_rows,
                                  std::array<Ipopt::Index, 5>{{0, 1, 2, 2, 2}});
          success *= checkEntries(hessian_columns,
                                  std::array<Ipopt::Index, 5>{{0, 1, 0, 1, 2}});

          const std::array<Ipopt::Number, 3> variables{{2.0, 3.0, 4.0}};
          const std::array<Ipopt::Number, 2> multipliers{{5.0, 7.0}};
          std::array<Ipopt::Number, 5>       jacobian_values{};
          std::array<Ipopt::Number, 5>       hessian_values{};
          std::array<Ipopt::Number, 3>       objective_gradient{};
          std::array<Ipopt::Number, 2>       constraint_values{};
          Ipopt::Number                      objective{0};

          success *= adapter.eval_f(
              variable_count, variables.data(), true, objective);
          success *= isEqual(objective,
                             static_cast<Ipopt::Number>(179.0 / 6.0),
                             tolerance());
          success *= adapter.eval_grad_f(variable_count,
                                         variables.data(),
                                         false,
                                         objective_gradient.data());
          success *= checkEntries(objective_gradient,
                                  std::array<Ipopt::Number, 3>{{4.0, 5.0, 21.0}});
          success *= adapter.eval_g(variable_count,
                                    variables.data(),
                                    false,
                                    constraint_count,
                                    constraint_values.data());
          success *= checkEntries(constraint_values,
                                  std::array<Ipopt::Number, 2>{{33.0, 7.0}});

          success *= adapter.eval_jac_g(variable_count,
                                        variables.data(),
                                        true,
                                        constraint_count,
                                        jacobian_nonzeros,
                                        nullptr,
                                        nullptr,
                                        jacobian_values.data());
          success *= checkEntries(jacobian_values,
                                  std::array<Ipopt::Number, 5>{{4.0, 4.0, 13.0, 3.0, 1.0}});

          success *= adapter.eval_h(variable_count,
                                    variables.data(),
                                    false,
                                    2.0,
                                    constraint_count,
                                    multipliers.data(),
                                    true,
                                    hessian_nonzeros,
                                    nullptr,
                                    nullptr,
                                    hessian_values.data());
          success *= checkEntries(hessian_values,
                                  std::array<Ipopt::Number, 5>{{2.0, 9.0, 7.0, 7.0, 28.0}});
          success *= !adapter.eval_h(variable_count,
                                     variables.data(),
                                     false,
                                     2.0,
                                     constraint_count,
                                     nullptr,
                                     false,
                                     hessian_nonzeros,
                                     nullptr,
                                     nullptr,
                                     hessian_values.data());

          MissingJacobianModel<ScalarT, IdxT> missing_jacobian;
          success *= missing_jacobian.allocate() == 0;
          success *= missing_jacobian.initialize() == 0;
          success *= GridKit::Testing::throws<std::invalid_argument>(
              [&missing_jacobian]()
              {
                AdapterT unavailable(missing_jacobian);
              });

          MissingHessianModel<ScalarT, IdxT> missing_hessian;
          success *= missing_hessian.allocate() == 0;
          success *= missing_hessian.initialize() == 0;
          success *= GridKit::Testing::throws<std::invalid_argument>(
              [&missing_hessian]()
              {
                AdapterT unavailable(missing_hessian);
              });

          return success.report(__func__);
        }

        TestOutcome emptyExactHessian()
        {
          TestStatus                 success = true;
          AffineModel<ScalarT, IdxT> model;
          success *= model.allocate() == 0;
          success *= model.initialize() == 0;

          AdapterT                    adapter(model);
          Ipopt::Index                variable_count{0};
          Ipopt::Index                constraint_count{0};
          Ipopt::Index                jacobian_nonzeros{0};
          Ipopt::Index                hessian_nonzeros{0};
          Ipopt::TNLP::IndexStyleEnum index_style{};
          success *= adapter.get_nlp_info(variable_count,
                                          constraint_count,
                                          jacobian_nonzeros,
                                          hessian_nonzeros,
                                          index_style);
          success *= variable_count == 1;
          success *= constraint_count == 0;
          success *= jacobian_nonzeros == 0;
          success *= hessian_nonzeros == 0;

          success *= adapter.eval_h(variable_count,
                                    nullptr,
                                    false,
                                    1.0,
                                    constraint_count,
                                    nullptr,
                                    false,
                                    hessian_nonzeros,
                                    nullptr,
                                    nullptr,
                                    nullptr);

          const std::array<Ipopt::Number, 1> variables{{0.25}};
          Ipopt::Number                      unused_value{0.0};
          success *= adapter.eval_h(variable_count,
                                    variables.data(),
                                    true,
                                    1.0,
                                    constraint_count,
                                    nullptr,
                                    false,
                                    hessian_nonzeros,
                                    nullptr,
                                    nullptr,
                                    &unused_value);
          success *= model.hessianEvaluationCount() == 1;
          return success.report(__func__);
        }

        TestOutcome secondOrderDerivativeTest()
        {
          TestStatus success = true;
          ModelT     model;
          success *= model.allocate() == 0;
          success *= model.initialize() == 0;

          Ipopt::SmartPtr<Ipopt::IpoptApplication> application  = IpoptApplicationFactory();
          success                                              *= application->Options()->SetStringValue("hessian_approximation", "exact");
          success                                              *= application->Options()->SetStringValue("derivative_test", "second-order");
          success                                              *= application->Options()->SetStringValue("derivative_test_print_all", "yes");
          success                                              *= application->Options()->SetNumericValue("derivative_test_perturbation", 1.0e-8);
          success                                              *= application->Options()->SetNumericValue("derivative_test_tol", 1.0e-6);
          success                                              *= application->Options()->SetNumericValue("tol", 1.0e-10);
          success                                              *= application->Options()->SetIntegerValue("print_level", 5);
          success                                              *= application->Options()->SetIntegerValue("max_iter", 100);
          success                                              *= application->Options()->SetIntegerValue("mumps_print_level", 0);

          auto status  = application->Initialize();
          success     *= status == Ipopt::Solve_Succeeded;
          if (status == Ipopt::Solve_Succeeded)
          {
            Ipopt::SmartPtr<Ipopt::TNLP> problem  = new AdapterT(model);
            status                                = application->OptimizeTNLP(problem);
            success                              *= status == Ipopt::Solve_Succeeded;
          }

          success *= checkVector(model.variables(),
                                 std::array<RealT, 3>{{1.0, 1.0, 1.0}},
                                 static_cast<RealT>(1.0e-7));
          success *= model.hessianEvaluationCount() > 0;
          success *= model.evaluateObjective() == 0;
          success *= isEqual(model.objective(),
                             static_cast<RealT>(-25.0 / 6.0),
                             static_cast<RealT>(1.0e-8));
          success *= model.evaluateConstraints() == 0;
          success *= checkVector(model.constraints(),
                                 std::array<RealT, 2>{{0.0, 0.0}},
                                 static_cast<RealT>(1.0e-8));
          return success.report(__func__);
        }

      private:
        static constexpr RealT tolerance()
        {
          return static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();
        }

        template <typename ValueT, std::size_t Count>
        static bool checkEntries(const std::array<ValueT, Count>& actual,
                                 const std::array<ValueT, Count>& expected)
        {
          for (std::size_t entry = 0; entry < Count; ++entry)
          {
            if constexpr (std::is_floating_point_v<ValueT>)
            {
              if (!isEqual(actual[entry], expected[entry], static_cast<ValueT>(tolerance())))
              {
                return false;
              }
            }
            else if (actual[entry] != expected[entry])
            {
              return false;
            }
          }
          return true;
        }

        template <std::size_t Count>
        static bool checkVector(const typename ModelT::VectorT& vector,
                                const std::array<RealT, Count>& expected,
                                RealT                           comparison_tolerance = tolerance())
        {
          if (vector.getSize() != static_cast<IdxT>(Count))
          {
            return false;
          }

          const auto* values = vector.getData(memory::HOST);
          for (std::size_t entry = 0; entry < Count; ++entry)
          {
            if (!isEqual(static_cast<RealT>(values[entry]),
                         expected[entry],
                         comparison_tolerance))
            {
              return false;
            }
          }
          return true;
        }
      };
    } // namespace Optimization
  } // namespace Testing
} // namespace GridKit
