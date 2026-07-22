#pragma once

#include <GridKit/Testing/Testing.hpp>

#include <zmq.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class ZMQTests
    {
    public:
      ZMQTests()  = default;
      ~ZMQTests() = default;

      TestOutcome basic()
      {
        TestStatus success = true;

        zmq::context_t ctx;
        zmq::socket_t  socket(ctx, zmq::socket_type::req);

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
