#pragma once

#include <sstream>
#include <string>
#include <vector>

#include <GridKit/Utilities/String.hpp>

namespace GridKit
{
  namespace Testing
  {
    /**
     * @brief Parse a line of tokens into a vector of values
     */
    template <typename T = std::string>
    class Tokenizer
    {
    public:
      Tokenizer() = delete;

      explicit Tokenizer(const std::string& in, char delimiter = ' ')
      {
        std::istringstream iss(in);
        for (std::string item; std::getline(iss, item, delimiter);)
        {
          if constexpr (std::is_same_v<T, std::string>)
          {
            tokens_.push_back(GridKit::Utilities::strip(item));
          }
          else
          {
            tokens_.push_back(GridKit::Utilities::parse<T>(item));
          }
        }
      }

      const std::vector<T>& operator()() const
      {
        return tokens_;
      }

    private:
      std::vector<T> tokens_;
    };

  } // namespace Testing
} // namespace GridKit
