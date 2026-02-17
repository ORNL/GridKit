#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <span>
#include <string>
#include <vector>

#include <GridKit/Testing/CSV.hpp>
#include <GridKit/Testing/OutputAtTime.hpp>
#include <GridKit/Testing/Tokenizer.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace Testing
  {

    using Log = ::GridKit::Utilities::Logger;

    std::ifstream checkOpenFile(const std::string& f)
    {
      if (!std::filesystem::exists(f))
      {
        Log::error() << "\n\tFile \"" << f << "\" does not exist." << std::endl;
        throw std::runtime_error("file not found");
      }
      std::ifstream ifs(f);
      if (!ifs)
      {
        Log::error() << "\n\tCould not open file \"" << f
                     << "\" for input." << std::endl;
        throw std::runtime_error("Failed to open file");
      }
      return ifs;
    }

    ErrorSet compareCSV(const std::string& f_a, const std::string& f_b)
    {
      // Open files and read labels
      auto ifs_a = checkOpenFile(f_a);
      auto ifs_b = checkOpenFile(f_b);

      std::string line_a;
      std::getline(ifs_a, line_a);
      auto labels_a = Tokenizer<>(line_a, ',')();

      std::string line_b;
      std::getline(ifs_b, line_b);
      auto labels_b = Tokenizer<>(line_b, ',')();

      if (labels_a.size() != labels_b.size())
      {
        Log::error() << "Files \"" << f_a << "\" and \"" << f_b
                     << "\" have different number of variables." << std::endl;
      }

      // Create error set
      auto err = ErrorSet(std::span{next(begin(labels_a)), end(labels_a)});

      while (ifs_a && ifs_b)
      {
        std::getline(ifs_a, line_a);
        std::getline(ifs_b, line_b);
        if (!(ifs_a && ifs_b))
        {
          if (ifs_a || ifs_b)
          {
            Log::error() << "Files \"" << f_a << "\" and \"" << f_b
                         << "\" have different lengths." << std::endl;
          }
          break;
        }

        OutputAtTime d_a(Tokenizer<double>(line_b, ',')());
        OutputAtTime d_b(Tokenizer<double>(line_a, ',')());

        err.push(d_a - d_b);
      }

      err.wrap();

      return err;
    }

  } // namespace Testing
} // namespace GridKit
