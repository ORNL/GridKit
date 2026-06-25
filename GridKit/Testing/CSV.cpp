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

    std::unique_ptr<ErrorSet> compareCSV(const std::string& f_test,
                                         const std::string& f_ref,
                                         ErrorType          error_type,
                                         double             abs_threshold)
    {
      // Open files and read labels
      auto ifs_test = checkOpenFile(f_test);
      auto ifs_ref  = checkOpenFile(f_ref);

      std::string line_test;
      std::getline(ifs_test, line_test);
      auto labels_a = Tokenizer<>(line_test, ',')();

      std::string line_ref;
      std::getline(ifs_ref, line_ref);
      auto labels_b = Tokenizer<>(line_ref, ',')();

      if (labels_a.size() != labels_b.size())
      {
        Log::error() << "Files \"" << f_test << "\" and \"" << f_ref
                     << "\" have different number of variables." << std::endl;
      }

      // Create error set
      auto var_labels = std::span{next(begin(labels_a)), end(labels_a)};
      auto err        = makeErrorSet(error_type, var_labels, abs_threshold);

      while (ifs_test && ifs_ref)
      {
        std::getline(ifs_test, line_test);
        std::getline(ifs_ref, line_ref);
        if (!(ifs_test && ifs_ref))
        {
          if (ifs_test || ifs_ref)
          {
            Log::error() << "Files \"" << f_test << "\" and \"" << f_ref
                         << "\" have different lengths." << std::endl;
          }
          break;
        }

        OutputAtTime d_test(Tokenizer<double>(line_test, ',')());
        OutputAtTime d_ref(Tokenizer<double>(line_ref, ',')());

        err->push(d_test, d_ref);
      }

      err->wrap();

      return err;
    }

  } // namespace Testing
} // namespace GridKit
