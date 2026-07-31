#include "SystemDataTests.hpp"

int main()
{
  GridKit::Testing::TestingResults  result;
  GridKit::Testing::SystemDataTests test;

  result += test.parseCompleteSystem();
  result += test.parseRepositoryCases();
  result += test.rejectInvalidDocuments();
  result += test.rejectMissingFile();

  return result.summary();
}
