gridkit_add_library(
  testing
  SOURCES Testing.cpp CSV.cpp
  HEADERS Testing.hpp
          TestHelpers.hpp
          Tokenizer.hpp
          OutputAtTime.hpp
          CSV.hpp
          AggregateErrors.hpp
          CommandLine.hpp
  LINK_LIBRARIES PUBLIC GridKit::Utilities)
