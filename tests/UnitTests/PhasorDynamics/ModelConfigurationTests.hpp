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
      negative,
      infinite,
      integer_flag,
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
      /// Omitted parameters keep the default, integers widen to reals, and
      /// invalid values are rejected naming the model and parameter.
      TestOutcome parameterReader()
      {
        TestStatus success = true;

        DataT data;
        data.parameters = {{Params::integer, static_cast<IdxT>(2)},
                           {Params::negative, -ONE<RealT>},
                           {Params::infinite, std::numeric_limits<RealT>::infinity()},
                           {Params::integer_flag, static_cast<IdxT>(1)}};

        const Model::ParameterReader reader(data, "Test");
        RealT                        real = ONE<RealT>;

        success *= !reader.loadReal(Params::omitted, real);
        success *= (real == ONE<RealT>);
        success *= (reader.requireReal(Params::integer) == TWO<RealT>);
        success *= throws<std::invalid_argument>([&]()
                                                 { reader.loadNonNegative(Params::negative, real); });
        success *= throws<std::invalid_argument>([&]()
                                                 { reader.loadReal(Params::infinite, real); });
        success *= throws<std::invalid_argument>([&]()
                                                 { reader.requireBool(Params::integer_flag); });

        try
        {
          reader.requireReal(Params::omitted);
          success = false;
        }
        catch (const std::invalid_argument& error)
        {
          success *= (std::string_view(error.what()) == "Test: parameter 'omitted' is required");
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
