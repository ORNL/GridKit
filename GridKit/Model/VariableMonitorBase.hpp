#pragma once

#include <ostream>
#include <string>

namespace GridKit
{
  namespace Model
  {
    enum class VariableMonitorFormat
    {
      CSV,
      JSON,
      YAML
    };

    class VariableMonitorBase
    {
    public:
      struct Csv
      {
        std::string delim{","};
      };

      struct Json
      {
        mutable bool after_first{false};
      };

      struct Yaml
      {
      };

      using Format = VariableMonitorFormat;

      struct SinkSpec
      {
        std::string file_name;
        Format      format;
        std::string delim;
      };

      virtual ~VariableMonitorBase()
      {
      }

      virtual bool empty() const = 0;

      virtual void printHeader(std::ostream&, Csv) const = 0;
      virtual void print(std::ostream&, Csv) const       = 0;

      virtual void printFooter(std::ostream&, Csv) const
      {
      }

      virtual void printHeader(std::ostream&, Json) const
      {
      }

      virtual void print(std::ostream&, Json) const = 0;

      virtual void printFooter(std::ostream&, Json) const
      {
      }

      virtual void printHeader(std::ostream&, Yaml) const
      {
      }

      virtual void print(std::ostream&, Yaml) const = 0;

      virtual void printFooter(std::ostream&, Yaml) const
      {
      }
    };

  } // namespace Model
} // namespace GridKit
