#pragma once

#include <array>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace GridKit
{
  namespace Utilities
  {

    class ArgValue
    {
    public:
      ArgValue() = default;

      template <typename T>
      ArgValue(T&& val)
        : value_((std::stringstream() << val).str())
      {
      }

      template <typename T>
      ArgValue& operator=(T&& val)
      {
        value_ = (std::stringstream() << val).str();
        return *this;
      }

      template <typename T>
      T as() const
      {
        T ret;
        std::stringstream(value_) >> ret;
        return ret;
      }

      const std::string& operator()() const
      {
        return value_;
      }

      const std::string& get() const
      {
        return value_;
      }

      bool empty() const
      {
        return value_.empty();
      }

    private:
      std::string value_{};
    };

    class ArgVector
    {
    public:
      ArgVector() = default;

      template <typename T>
      ArgVector(T&& val)
        : vec{val}
      {
      }

      template <typename T, std::size_t N = 1>
      std::array<T, N> as() const
      {
        assert(vec.size() == N);
        std::array<T, N> ret;
        for (std::size_t i = 0; i < N; ++i)
        {
          ret[i] = vec[i].as<T>();
        }
        return ret;
      }

      const ArgValue& operator()() const
      {
        return vec[0];
      }

      const ArgValue& operator[](std::size_t i) const
      {
        return vec[i];
      }

      bool empty() const
      {
        return vec.empty();
      }

      std::size_t size() const
      {
        return vec.size();
      }

    private:
      std::vector<ArgValue> vec;

      friend class CliArgs;
    };

    std::ostream& operator<<(std::ostream& os, const ArgVector& av);

    struct Argument
    {
      enum class Type
      {
        String,
        Real,
        Integer,
        Boolean,
        Unspecified
      };

      std::vector<std::string> name;
      std::string              help{};
      bool                     required{false};
      Type                     type{Type::Unspecified};
      bool                     flag{false};
      ArgVector                defaults{};
      std::size_t              nargs{1};
      ArgVector                values{};
    };

    std::ostream& operator<<(std::ostream& os, const Argument& arg);

    class CliArgs
    {
    public:
      CliArgs(std::initializer_list<Argument> args);

      void parseArgs(int argc, const char* argv[]);

      void printHelp(std::ostream& os = std::cout) const;

      const ArgVector& operator[](const std::string& name) const
      {
        auto arg = name_map_.at(name);
        if (arg->values.empty() && !arg->defaults.empty())
        {
          return arg->defaults;
        }
        return arg->values;
      }

    private:
      void      mapName(const std::string& name, Argument* arg);
      void      setupArg(Argument& arg);
      Argument* findArg(const std::string& raw);

      std::string                                app_name_;
      std::vector<Argument>                      table_;
      std::unordered_map<std::string, Argument*> name_map_;

      friend std::ostream& operator<<(std::ostream& os, const CliArgs& args);
    };

  } // namespace Utilities
} // namespace GridKit
