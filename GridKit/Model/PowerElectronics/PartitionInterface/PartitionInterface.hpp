

#pragma once

#include <algorithm>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

namespace GridKit
{
  /*!
   * @brief Base class for partition interface components.
   *
   * PartitionInterface provides the infrastructure needed by components
   * that exchange information across partition boundaries during
   * co-simulation. It stores the values of external variables received
   * from neighboring partitions together with the corresponding global
   * indices in the original unpartitioned system.
   *
   * Derived classes are responsible for implementing the specific
   * interface behavior associated with buses, circuit components, or
   * other partition boundary entities.
   */
  template <class ScalarT, typename IdxT>
  class PartitionInterface : public CircuitComponent<ScalarT, IdxT>
  {
  public:
    using RealT      = typename Model::Evaluator<ScalarT, IdxT>::RealT;
    using MatrixT    = typename Model::Evaluator<ScalarT, IdxT>::MatrixT;
    using CsrMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CsrMatrixT;

    PartitionInterface() = default;

    ~PartitionInterface() = default;

    int setExternalDataY(const std::vector<RealT>& data)
    {
      std::copy(data.begin(), data.end(), external_data_y_.begin());

      return 0;
    }

    int setExternalDataYP(const std::vector<RealT>& data)
    {
      std::copy(data.begin(), data.end(), external_data_yp_.begin());

      return 0;
    }

    std::vector<RealT>& getExternalDataY()
    {
      return external_data_y_;
    }

    std::vector<RealT>& getExternalDataYP()
    {
      return external_data_yp_;
    }

    std::vector<IdxT>& getPartitionExternalIndices()
    {
      return interface_partition_externals_;
    }

  protected:
    std::vector<RealT> external_data_y_;
    std::vector<RealT> external_data_yp_;
    std::vector<IdxT>  interface_partition_externals_;
  };

} // namespace GridKit
