#include "BusTests.hpp"

int main()
{
    using namespace GridKit;
    using namespace GridKit::Testing;

    GridKit::Testing::TestingResults result; 
    GridKit::Testing::BusTests test;

    result += test.smoke();

    return result.summary();
}