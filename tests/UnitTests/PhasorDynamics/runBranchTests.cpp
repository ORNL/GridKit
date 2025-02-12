#include "BranchTests.hpp"

int main()
{
    using namespace GridKit;
    using namespace GridKit::Testing;

    GridKit::Testing::TestingResults result; 
    GridKit::Testing::BranchTests<double, size_t> test;

    result += test.smoke();

    return result.summary();
}