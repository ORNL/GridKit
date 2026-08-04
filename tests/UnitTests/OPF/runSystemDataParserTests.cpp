#include "SystemDataParserTests.hpp"

int main()
{
  GridKit::Testing::TestingResults        result;
  GridKit::Testing::SystemDataParserTests test;

  result += test.parseCompleteSystem();
  result += test.parseOptionalFieldsAndDefaults();
  result += test.rejectUnknownFields();
  result += test.rejectInvalidTypesAndVersions();
  result += test.rejectMissingRequiredFields();
  result += test.rejectMissingFile();
  result += test.parseCuratedFixtures();

  return result.summary();
}
