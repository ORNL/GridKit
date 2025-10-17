#pragma once

#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

namespace GridKit
{
  namespace Utilities
  {

    /**
     * @brief Parser for command line input
     *
     * @note Based on StackOverflow answer by Luca Davanzo
     */
    class CliOptions
    {
    public:
      CliOptions(int argc, const char* argv[]);
      virtual ~CliOptions();
      std::string        getAppName() const;
      bool               hasKey(const std::string&) const;
      const std::string& getOptionFromKey(const std::string&) const;
      void               printOptionsList(std::ostream& os = std::cout) const;

      template <typename T = std::string>
      T get(const std::string& key) const
      {
        T ret;
        std::stringstream(this->getOptionFromKey(key)) >> ret;
        return ret;
      }

    private:
      using OptionsList = std::map<std::string, std::string>;
      void               parse();
      const char* const* begin() const;
      const char* const* end() const;
      const char* const* last() const;
      OptionsList        options_;
      int                argc_;
      const char**       argv_;
      std::string        app_name_;
    };

  } // namespace Utilities
} // namespace GridKit
