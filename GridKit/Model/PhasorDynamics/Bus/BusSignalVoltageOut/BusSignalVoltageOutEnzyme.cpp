/**
 * @file BusSignalVoltageOutEnzyme.cpp
 * @author Slaven Peles (peless@ornl.gov)
 */

#include "BusSignalVoltageOutImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation for the Enzyme build.
     *
     * The first four entries are zero placeholders for the bus voltage
     * columns, as in Bus. They provide the indices for entries that other
     * components contribute to and that are later deduplicated.
     *
     * One additional entry per connected signal inlet holds the derivative of
     * the residual with respect to the signal variable, which is one. When
     * the signal has no valid system variable index, the entry is stored as
     * a zero duplicate of the bus-voltage placeholder so that the sparsity
     * pattern stays fixed.
     *
     * @return int - error code
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageOut<scalar_type, index_type>::evaluateJacobian()
    {
      constexpr IdxT num_bus_entries = 4;

      const auto& ir_port = ports_.in.template port<BusSignalVoltageOutInputs::ir>();
      const auto& ii_port = ports_.in.template port<BusSignalVoltageOutInputs::ii>();

      if (J_rows_buffer_ == nullptr)
      {
        nnz_ = num_bus_entries;
        if (ir_port.connected())
        {
          ++nnz_;
        }
        if (ii_port.connected())
        {
          ++nnz_;
        }

        const auto num_entries = static_cast<size_t>(nnz_);
        J_rows_buffer_         = new IdxT[num_entries];
        J_cols_buffer_         = new IdxT[num_entries];
        J_vals_buffer_         = new RealT[num_entries];
      }

      J_rows_buffer_[0] = residual_indices_.at(0);
      J_rows_buffer_[1] = residual_indices_.at(0);
      J_rows_buffer_[2] = residual_indices_.at(1);
      J_rows_buffer_[3] = residual_indices_.at(1);
      J_cols_buffer_[0] = variable_indices_.at(0);
      J_cols_buffer_[1] = variable_indices_.at(1);
      J_cols_buffer_[2] = variable_indices_.at(0);
      J_cols_buffer_[3] = variable_indices_.at(1);
      J_vals_buffer_[0] = 0.0;
      J_vals_buffer_[1] = 0.0;
      J_vals_buffer_[2] = 0.0;
      J_vals_buffer_[3] = 0.0;

      IdxT k = num_bus_entries;

      auto add_signal_entry = [&](const auto& port, size_t residual)
      {
        if (!port.connected())
        {
          return;
        }
        const IdxT signal_index = port.linked() ? port.signalVariableIndex() : INVALID_INDEX<IdxT>;
        J_rows_buffer_[k]       = residual_indices_.at(residual);
        if (signal_index != INVALID_INDEX<IdxT>)
        {
          J_cols_buffer_[k] = signal_index;
          J_vals_buffer_[k] = 1.0;
        }
        else
        {
          J_cols_buffer_[k] = variable_indices_.at(residual);
          J_vals_buffer_[k] = 0.0;
        }
        ++k;
      };

      add_signal_entry(ir_port, 0);
      add_signal_entry(ii_port, 1);

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class BusSignalVoltageOut<double, long int>;
    template class BusSignalVoltageOut<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
