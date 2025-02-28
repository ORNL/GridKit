#include "SynchronousMachineTests.hpp"

int main()
{
    using namespace GridKit;
    using namespace GridKit::Testing;

    GridKit::Testing::TestingResults                          result;
    GridKit::Testing::SynchronousMachineTests<double, size_t> test;

    result += test.constructor();
    result += test.accessors();
    result += test.residual();

    return result.summary();
}
