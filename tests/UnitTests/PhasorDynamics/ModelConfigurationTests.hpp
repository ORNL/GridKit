#pragma once

#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace Testing
  {
    /// One parameter per ParameterReader case, named for its value.
    enum class ReaderTestParameters : std::size_t
    {
      omitted,
      integer,
      real,
      zero,
      negative,
      infinite,
      not_a_number,
      integer_flag,
      false_flag,
      true_flag,
      text,
    };

    template <typename RealT, typename IdxT>
    class ModelConfigurationTests
    {
    private:
      using Params = ReaderTestParameters;
      using DataT  = PhasorDynamics::ComponentData<RealT,
                                                   IdxT,
                                                   Params,
                                                   Utilities::EmptyEnum,
                                                   Utilities::EmptyEnum,
                                                   Utilities::EmptyEnum,
                                                   Utilities::EmptyEnum>;

    public:
      /// Omitted parameters return no value, integers widen to reals, and
      /// invalid values are rejected naming the model and parameter.
      TestOutcome parameterReader()
      {
        TestStatus success = true;

        DataT data;
        data.parameters = {{Params::integer, static_cast<IdxT>(2)},
                           {Params::real, static_cast<RealT>(0.5)},
                           {Params::zero, ZERO<RealT>},
                           {Params::negative, -ONE<RealT>},
                           {Params::infinite, std::numeric_limits<RealT>::infinity()},
                           {Params::not_a_number, std::numeric_limits<RealT>::quiet_NaN()},
                           {Params::integer_flag, static_cast<IdxT>(1)},
                           {Params::false_flag, false},
                           {Params::true_flag, true},
                           {Params::text, std::string("2")}};

        const Model::ParameterReader reader(data, "Test");

        success *= !reader.readReal(Params::omitted).has_value();
        success *= !reader.readNonNegative(Params::omitted).has_value();
        success *= !reader.readBool(Params::omitted).has_value();
        success *= !reader.readInteger(Params::omitted).has_value();
        success *= !reader.readString(Params::omitted).has_value();
        success *= (reader.readReal(Params::omitted).value_or(ONE<RealT>) == ONE<RealT>);
        success *= reader.readBool(Params::omitted).value_or(true);

        success *= (reader.readReal(Params::integer) == TWO<RealT>);
        success *= (reader.readReal(Params::real) == static_cast<RealT>(0.5));
        success *= (reader.readReal(Params::negative) == -ONE<RealT>);
        success *= (reader.readNonNegative(Params::zero) == ZERO<RealT>);
        success *= (reader.readNonNegative(Params::integer) == TWO<RealT>);
        success *= (reader.readInteger(Params::integer) == static_cast<IdxT>(2));
        success *= (reader.readBool(Params::false_flag) == false);
        success *= (reader.readBool(Params::true_flag) == true);
        success *= !reader.readBool(Params::false_flag).value_or(true);
        success *= (reader.readString(Params::text) == "2");

        success *= throws<std::invalid_argument>([&]()
                                                 { static_cast<void>(reader.readReal(Params::infinite)); });
        success *= throws<std::invalid_argument>([&]()
                                                 { static_cast<void>(reader.readReal(Params::not_a_number)); });
        success *= throws<std::invalid_argument>([&]()
                                                 { static_cast<void>(reader.readReal(Params::text)); });
        success *= throws<std::invalid_argument>([&]()
                                                 { static_cast<void>(reader.readReal(Params::true_flag)); });
        success *= throws<std::invalid_argument>([&]()
                                                 { static_cast<void>(reader.readInteger(Params::real)); });
        success *= throws<std::invalid_argument>([&]()
                                                 { static_cast<void>(reader.readInteger(Params::zero)); });
        success *= throws<std::invalid_argument>([&]()
                                                 { static_cast<void>(reader.readInteger(Params::true_flag)); });
        success *= throws<std::invalid_argument>([&]()
                                                 { static_cast<void>(reader.readNonNegative(Params::negative)); });
        success *= throws<std::invalid_argument>([&]()
                                                 { static_cast<void>(reader.readString(Params::integer)); });

        try
        {
          static_cast<void>(reader.readBool(Params::integer_flag));
          success = false;
        }
        catch (const std::invalid_argument& error)
        {
          success *= (std::string_view(error.what()) == "Test: parameter 'integer_flag' must be boolean");
        }

        return success.report(__func__);
      }

      /// Checks record only failed conditions, prefixed by the model name.
      TestOutcome configurationChecks()
      {
        TestStatus success = true;

        Model::ConfigurationChecks checks("Test");
        checks.check(true, "never recorded");
        success *= checks.passed();

        checks.check(false, "x must be positive");
        success *= !checks.passed();
        success *= (checks.errors() == std::vector<std::string>{"Test: x must be positive"});

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
