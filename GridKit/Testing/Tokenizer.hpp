#pragma once

#include <sstream>
#include <string>
#include <vector>

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
          std::istringstream(item) >> tokens_.emplace_back();
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
