#include "SparsityPatternTests.hpp"

int main()
{
    using namespace GridKit;
    using namespace GridKit::Testing;

    GridKit::Testing::TestingResults result; 
    GridKit::Testing::SparsityPatternTests<double, size_t> test;

    result += test.variable();

    return result.summary();
}
