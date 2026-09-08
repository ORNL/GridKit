#pragma once

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    class Delay;

    /**
     * @brief Accepted-history discontinuities that restart the integrator.
     *
     * Delay is the only owner. Every other EMT model is smooth by design.
     * The private constructor makes any other derivation a compile error.
     */
    template <typename RealT>
    class HistoryDiscontinuity
    {
    public:
      virtual RealT nextHistoryDiscontinuity(RealT after) const = 0;

    protected:
      ~HistoryDiscontinuity() = default;

    private:
      HistoryDiscontinuity() = default;

      template <typename scalar_type, typename index_type>
      friend class Delay;
    };
  } // namespace EMT
} // namespace GridKit
