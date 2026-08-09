/**
 * @file SignalNodeJunctionTests.hpp
 * @brief Unit tests for algebraic signal-node junctions.
 */
#pragma once

#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeJunction.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class SignalNodeJunctionTests
    {
    public:
      using ScalarT       = scalar_type;
      using IdxT          = index_type;
      using SignalT       = PhasorDynamics::SignalNode<ScalarT, IdxT>;
      using JunctionT     = PhasorDynamics::SignalNodeJunction<ScalarT, IdxT>;
      using JunctionInput = typename SignalT::JunctionInput;
      using RealT         = typename JunctionT::RealT;

      static constexpr RealT kTol = std::numeric_limits<RealT>::epsilon();

      /**
       * @brief Check forward initialization, the algebraic residual, and its Jacobian.
       */
      TestOutcome algebraicEquation()
      {
        TestStatus success = true;

        const RealT bias{0.5};
        const RealT gain0{2.0};
        const RealT gain1{-0.5};
        const RealT prescribed_gain{0.25};

        ScalarT input0_value{1.5};
        ScalarT input1_value{-2.0};
        ScalarT prescribed_value{4.0};

        IdxT input0_index{4};
        IdxT input1_index{6};
        IdxT prescribed_index{INVALID_INDEX<IdxT>};

        SignalT input0;
        SignalT input1;
        SignalT prescribed;
        SignalT output;

        input0.set(&input0_value, &input0_index);
        input1.set(&input1_value, &input1_index);
        prescribed.set(&prescribed_value, &prescribed_index);

        std::vector<JunctionInput> inputs{
            {&input0, gain0},
            {&input1, gain1},
            {&prescribed, prescribed_gain}};
        JunctionT junction(&output, bias, 0, std::move(inputs));

        success *= junction.allocate() == 0;

        const IdxT output_index{8};
        const IdxT residual_index{9};
        junction.setVariableIndex(0, output_index);
        junction.setResidualIndex(0, residual_index);

        success *= junction.verify() == 0;
        success *= junction.initialize() == 0;

        const ScalarT expected = static_cast<ScalarT>(bias)
                                 + gain0 * input0_value
                                 + gain1 * input1_value
                                 + prescribed_gain * prescribed_value;
        success *= isEqual(output.read(), expected, static_cast<ScalarT>(kTol));
        success *= isEqual(output.junctionValue(), expected, static_cast<ScalarT>(kTol));

        const ScalarT residual_offset{0.25};
        junction.y().getData()[0]  = expected + residual_offset;
        success                   *= junction.evaluateResidual() == 0;
        success                   *= isEqual(junction.getResidual().getData()[0],
                           residual_offset,
                           static_cast<ScalarT>(kTol));

        success        *= junction.evaluateJacobian() == 0;
        auto* jacobian  = junction.getCooJacobian();
        success        *= jacobian != nullptr;
        if (jacobian != nullptr)
        {
          success *= jacobian->getNnz() == 3;

          const auto* rows    = jacobian->getRowData();
          const auto* columns = jacobian->getColData();
          const auto* values  = jacobian->getValues();

          std::map<IdxT, RealT> jacobian_row;
          for (IdxT entry = 0; entry < jacobian->getNnz(); ++entry)
          {
            success                      *= rows[entry] == residual_index;
            jacobian_row[columns[entry]] += values[entry];
          }

          const std::map<IdxT, RealT> expected_row{
              {output_index, ONE<RealT>},
              {input0_index, -gain0},
              {input1_index, -gain1}};
          success *= isEqual(jacobian_row, expected_row, kTol);
        }

        return success.report(__func__);
      }

      /**
       * @brief Check that output initialization back-solves the designated input.
       */
      TestOutcome designatedInputInitialization()
      {
        TestStatus success = true;

        const RealT bias{0.5};
        const RealT fixed_gain{3.0};
        const RealT initialization_gain{-2.0};

        ScalarT fixed_value{2.0};
        ScalarT initialization_value{0.0};
        IdxT    fixed_index{1};
        IdxT    initialization_index{2};

        SignalT fixed_input;
        SignalT initialization_input;
        SignalT output;
        fixed_input.set(&fixed_value, &fixed_index);
        initialization_input.set(&initialization_value, &initialization_index);

        std::vector<JunctionInput> inputs{
            {&fixed_input, fixed_gain},
            {&initialization_input, initialization_gain}};
        JunctionT junction(&output, bias, 1, std::move(inputs));
        success *= junction.allocate() == 0;

        const ScalarT requested_output{5.5};
        const ScalarT expected_initialization_value =
            (requested_output - static_cast<ScalarT>(bias) - fixed_gain * fixed_value)
            / initialization_gain;

        output.init(requested_output);

        success *= isEqual(output.read(), requested_output, static_cast<ScalarT>(kTol));
        success *= isEqual(fixed_input.read(), fixed_value, static_cast<ScalarT>(kTol));
        success *= isEqual(initialization_input.read(),
                           expected_initialization_value,
                           static_cast<ScalarT>(kTol));
        success *= isEqual(output.junctionValue(), requested_output, static_cast<ScalarT>(kTol));

        return success.report(__func__);
      }

      /**
       * @brief Check local junction topology validation.
       */
      TestOutcome invalidConfiguration()
      {
        TestStatus success = true;

        SignalT input;

        success *= throws<std::invalid_argument>([&]()
                                                 {
                                                   std::vector<JunctionInput> inputs{{&input, ONE<RealT>}};
                                                   JunctionT junction(nullptr, ZERO<RealT>, 0, std::move(inputs)); });

        success *= throws<std::invalid_argument>([&]()
                                                 {
                                                   SignalT                   output;
                                                   std::vector<JunctionInput> inputs;
                                                   JunctionT junction(&output, ZERO<RealT>, 0, std::move(inputs)); });

        success *= throws<std::invalid_argument>([&]()
                                                 {
                                                   SignalT output;
                                                   std::vector<JunctionInput> inputs{{&input, ONE<RealT>}};
                                                   JunctionT junction(&output, ZERO<RealT>, 1, std::move(inputs)); });

        success *= throws<std::invalid_argument>([&]()
                                                 {
                                                   SignalT output;
                                                   std::vector<JunctionInput> inputs{{nullptr, ONE<RealT>}};
                                                   JunctionT junction(&output, ZERO<RealT>, 0, std::move(inputs)); });

        success *= throws<std::invalid_argument>([&]()
                                                 {
                                                   SignalT output;
                                                   std::vector<JunctionInput> inputs{{&output, ONE<RealT>}};
                                                   JunctionT junction(&output, ZERO<RealT>, 0, std::move(inputs)); });

        success *= throws<std::invalid_argument>([&]()
                                                 {
                                                   SignalT output;
                                                   std::vector<JunctionInput> inputs{
                                                       {&input, ONE<RealT>},
                                                       {&input, ONE<RealT>}};
                                                   JunctionT junction(&output, ZERO<RealT>, 0, std::move(inputs)); });

        success *= throws<std::invalid_argument>([&]()
                                                 {
                                                   SignalT output;
                                                   std::vector<JunctionInput> inputs{{&input, ZERO<RealT>}};
                                                   JunctionT junction(&output, ZERO<RealT>, 0, std::move(inputs)); });

        success *= throws<std::invalid_argument>([&]()
                                                 {
                                                   SignalT output;
                                                   std::vector<JunctionInput> inputs{{&input, ONE<RealT>}};
                                                   JunctionT junction(
                                                       &output,
                                                       std::numeric_limits<RealT>::quiet_NaN(),
                                                       0,
                                                       std::move(inputs)); });

        success *= throws<std::invalid_argument>([&]()
                                                 {
                                                   SignalT output;
                                                   std::vector<JunctionInput> inputs{
                                                       {&input, std::numeric_limits<RealT>::infinity()}};
                                                   JunctionT junction(&output, ZERO<RealT>, 0, std::move(inputs)); });

        ScalarT produced_value{ZERO<ScalarT>};
        IdxT    produced_index{0};
        SignalT produced_output;
        produced_output.set(&produced_value, &produced_index);
        std::vector<JunctionInput> produced_inputs{{&input, ONE<RealT>}};
        JunctionT                  produced_junction(
            &produced_output, ZERO<RealT>, 0, std::move(produced_inputs));
        success *= produced_junction.allocate() != 0;

        SignalT                    valid_output;
        std::vector<JunctionInput> valid_inputs{{&input, ONE<RealT>}};
        JunctionT                  valid_junction(&valid_output, ZERO<RealT>, 0, std::move(valid_inputs));
        success *= valid_output.isJunction();
        success *= throws<std::logic_error>([&]()
                                            {
                                              std::vector<JunctionInput> replacement{{&input, ONE<RealT>}};
                                              valid_output.configureJunction(ZERO<RealT>, 0, std::move(replacement)); });

        return success.report(__func__);
      }

      /**
       * @brief Check the runtime guard for an indirect initialization cycle.
       */
      TestOutcome initializationCycle()
      {
        TestStatus success = true;

        SignalT first_output;
        SignalT second_output;

        std::vector<JunctionInput> first_inputs{{&second_output, ONE<RealT>}};
        std::vector<JunctionInput> second_inputs{{&first_output, ONE<RealT>}};
        JunctionT                  first(&first_output, ZERO<RealT>, 0, std::move(first_inputs));
        JunctionT                  second(&second_output, ZERO<RealT>, 0, std::move(second_inputs));

        success *= first.allocate() == 0;
        success *= second.allocate() == 0;

        success *= throws<std::logic_error>([&]()
                                            { first_output.init(static_cast<ScalarT>(1.0)); });
        success *= throws<std::logic_error>([&]()
                                            { second_output.init(static_cast<ScalarT>(2.0)); });

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
