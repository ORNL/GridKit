#include <algorithm>
#include <cassert>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>

#include <GridKit/Utilities/CliArgs/CliArgs.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "CliArgsImpl.hpp"

namespace GridKit
{
  namespace Utilities
  {
    CliArgs::CliArgs(std::initializer_list<Option> opts)
      : pImpl_(std::make_unique<CliArgsImpl>(opts))
    {
    }

    CliArgs::~CliArgs()
    {
    }

    std::ostream& operator<<(std::ostream& os, const CliArgs& args)
    {
      os << *(args.pImpl_);
      return os;
    }

    void CliArgs::parseArgs(int argc, const char* argv[])
    {
      pImpl_->parseArgs(argc, argv);
    }

    void CliArgs::printUsage(std::ostream& os) const
    {
      pImpl_->printUsage(os);
    }

    void CliArgs::printHelp(std::ostream& os) const
    {
      pImpl_->printHelp(os);
    }

    const std::string& CliArgs::getAppName() const
    {
      return pImpl_->app_name_;
    }

    const ArgVector& CliArgs::operator[](const std::string& name) const
    {
      return (*pImpl_)[name];
    }

  } // namespace Utilities
} // namespace GridKit
