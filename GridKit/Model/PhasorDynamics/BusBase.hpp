#pragma once

#include <cassert>
#include <memory>
#include <set>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Constants.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>
#include <GridKit/Utilities/Errors.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /*!
     * @brief BusBase model implementation base class.
     *
     */
    template <typename scalar_type, typename index_type>
    class BusBase : public Model::Evaluator<scalar_type, index_type>
    {
    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using CsrMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CsrMatrixT;
      using CooMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CooMatrixT;
      using VectorT    = typename Model::Evaluator<ScalarT, IdxT>::VectorT;
      using BusTypeT   = typename BusData<RealT, IdxT>::BusType;
      using MonitorT   = Model::VariableMonitor<BusBase, BusData>;

      BusBase() = default;

      virtual ~BusBase();

      virtual int verify() const
      {
        return 0;
      }

      IdxT size() override final
      {
        return size_;
      }

      IdxT nnz() override final
      {
        return nnz_;
      }

      VectorT& y() override
      {
        return y_;
      }

      const VectorT& y() const override
      {
        return y_;
      }

      VectorT& yp() override
      {
        return yp_;
      }

      const VectorT& yp() const override
      {
        return yp_;
      }

      std::vector<bool>& tag() override
      {
        return tag_;
      }

      const std::vector<bool>& tag() const override
      {
        return tag_;
      }

      VectorT& absoluteTolerance() override
      {
        return abs_tol_;
      }

      const VectorT& absoluteTolerance() const override
      {
        return abs_tol_;
      }

      VectorT& getResidual() override
      {
        return f_;
      }

      const VectorT& getResidual() const override
      {
        return f_;
      }

      /**
       * @brief Bind this bus's state and residual vectors to the slice
       * [offset, offset + size()) of the system vectors.
       *
       * After binding, the bus reads and writes system storage directly
       * and allocate() will not allocate local vector data. Rebinding is
       * allowed and refreshes the aliases, e.g. after system storage is
       * reallocated when the topology changes.
       *
       * Only HOST data is bound because PhasorDynamics currently evaluates
       * models on the CPU. Supporting DEVICE data would also require sharing
       * the matching DEVICE pointer and keeping the HOST and DEVICE copies in
       * sync. This bind operation does neither, so DEVICE access is unsupported.
       *
       * @param[in] y       - System state vector
       * @param[in] yp      - System state derivative vector
       * @param[in] f       - System residual vector
       * @param[in] abs_tol - System absolute tolerance vector
       * @param[in] offset  - Position of this bus's slice in the system vectors
       *
       * @pre System vectors hold current HOST data of at least offset + size()
       * elements. This bus's vectors are unallocated or already bound.
       * @post allocated_ is true; y_, yp_, f_, and abs_tol_ alias system
       * storage; terminal accessors reference the newly bound storage.
       *
       * @return 0 if successful, non-zero otherwise.
       */
      int bind(VectorT& y, VectorT& yp, VectorT& f, VectorT& abs_tol, IdxT offset)
      {
        if (y.getSize() < offset + size_
            || yp.getSize() < offset + size_
            || f.getSize() < offset + size_
            || abs_tol.getSize() < offset + size_)
        {
          Log::error() << "BusBase::bind - system vectors are smaller than "
                       << "offset + size = " << offset + size_ << "\n";
          return 1;
        }

        auto* y_data       = y.getData(memory::HOST);
        auto* yp_data      = yp.getData(memory::HOST);
        auto* f_data       = f.getData(memory::HOST);
        auto* abs_tol_data = abs_tol.getData(memory::HOST);

        if (y_data == nullptr || yp_data == nullptr
            || f_data == nullptr || abs_tol_data == nullptr)
        {
          Log::error() << "BusBase::bind - system vector data is null or stale\n";
          return 1;
        }

        const int y_status       = y_.setData(y_data + offset, size_, memory::HOST);
        const int yp_status      = yp_.setData(yp_data + offset, size_, memory::HOST);
        const int f_status       = f_.setData(f_data + offset, size_, memory::HOST);
        const int abs_tol_status = abs_tol_.setData(abs_tol_data + offset, size_, memory::HOST);

        if (y_status != 0 || yp_status != 0 || f_status != 0 || abs_tol_status != 0)
        {
          Log::error() << "BusBase::bind - failed to bind vectors to system storage\n";
          return 1;
        }

        if (refreshTerminals() != 0)
        {
          Log::error() << "BusBase::bind - failed to refresh terminal storage\n";
          return 1;
        }

        allocated_ = true;
        return 0;
      }

      int setVariableIndex(IdxT local_index, IdxT global_index)
      {
        variable_indices_[static_cast<size_t>(local_index)] = global_index;
        return 0;
      }

      IdxT& getVariableIndex(IdxT local_index)
      {
        return variable_indices_[static_cast<size_t>(local_index)];
      }

      const std::vector<IdxT>& getVariableIndices() const
      {
        return variable_indices_;
      }

      int setResidualIndex(IdxT local_index, IdxT global_index)
      {
        residual_indices_[static_cast<size_t>(local_index)] = global_index;
        return 0;
      }

      IdxT& getResidualIndex(IdxT local_index)
      {
        return residual_indices_[static_cast<size_t>(local_index)];
      }

      const std::vector<IdxT>& getResidualIndices() const
      {
        return residual_indices_;
      }

      /// Pure virtual function, returns bus type (DEFAULT or SLACK).
      virtual BusTypeT BusType() const
      {
        return BusTypeT::DEFAULT;
      }

      CsrMatrixT* getCsrJacobian() const override
      {
        return csr_jac_;
      }

      CooMatrixT* getCooJacobian() const
      {
        return coo_jac_;
      }

      bool hasJacobian() override
      {
        return true;
      }

      void updateTime(RealT /* t */, RealT /* a */) override
      {
        // No time to update in bus models
      }

      /**
       * @pre Terminal storage has been established by Bus::allocate(),
       * BusBase::bind(), or the BusInfinite constructor.
       */
      ScalarT& Vr()
      {
        assert(Vr_ptr_ != nullptr);
        return *Vr_ptr_;
      }

      const ScalarT& Vr() const
      {
        assert(Vr_ptr_ != nullptr);
        return *Vr_ptr_;
      }

      ScalarT& Vi()
      {
        assert(Vi_ptr_ != nullptr);
        return *Vi_ptr_;
      }

      const ScalarT& Vi() const
      {
        assert(Vi_ptr_ != nullptr);
        return *Vi_ptr_;
      }

      ScalarT& Ir()
      {
        assert(Ir_ptr_ != nullptr);
        return *Ir_ptr_;
      }

      const ScalarT& Ir() const
      {
        assert(Ir_ptr_ != nullptr);
        return *Ir_ptr_;
      }

      ScalarT& Ii()
      {
        assert(Ii_ptr_ != nullptr);
        return *Ii_ptr_;
      }

      const ScalarT& Ii() const
      {
        assert(Ii_ptr_ != nullptr);
        return *Ii_ptr_;
      }

      virtual int setBusID(IdxT) = 0;

      virtual const IdxT busID() const
      {
        return bus_id_;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    protected:
      /**
       * @brief Refresh cached HOST terminal pointers after allocation or rebinding.
       */
      virtual int refreshTerminals() = 0;

      void setTerminals(ScalarT* Vr, ScalarT* Vi, ScalarT* Ir, ScalarT* Ii)
      {
        Vr_ptr_ = Vr;
        Vi_ptr_ = Vi;
        Ir_ptr_ = Ir;
        Ii_ptr_ = Ii;
      }

      /**
       * @brief Allocate this bus's state and residual vectors.
       */
      void allocateVectors(IdxT n)
      {

        y_.resize(n);
        yp_.resize(n);
        f_.resize(n);
        abs_tol_.resize(n);
      }

      IdxT              size_{0};
      IdxT              nnz_{0};
      /// Global (system-level) variable indices
      std::vector<IdxT> variable_indices_;
      /// Global (system-level) residual indices
      std::vector<IdxT> residual_indices_;

      VectorT           y_;
      VectorT           yp_;
      std::vector<bool> tag_;
      VectorT           abs_tol_;
      VectorT           f_;
      bool              allocated_{false};

      ScalarT* Vr_ptr_{nullptr};
      ScalarT* Vi_ptr_{nullptr};
      ScalarT* Ir_ptr_{nullptr};
      ScalarT* Ii_ptr_{nullptr};

      std::vector<ScalarT> g_;

      CsrMatrixT* csr_jac_{nullptr};
      CooMatrixT* coo_jac_{nullptr};

      //
      // Adjoint sensitivity members
      //

      std::vector<ScalarT> yB_{};
      std::vector<ScalarT> ypB_{};
      std::vector<ScalarT> fB_{};
      std::vector<ScalarT> gB_{};

      std::vector<ScalarT> param_{};
      std::vector<ScalarT> param_up_{};
      std::vector<ScalarT> param_lo_{};

      IdxT bus_id_{INVALID_INDEX<IdxT>};

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;

      using NotImplementedError = GridKit::Utilities::NotImplementedError;

    public:
      // TODO: evaluate how this complies with xSDK guidelines

      [[noreturn]] IdxT sizeQuadrature() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] IdxT sizeParams() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& yB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& yB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& ypB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& ypB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& param() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& param() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& param_up() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& param_up() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& param_lo() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& param_lo() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int initializeAdjoint() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateAdjointResidual() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateAdjointIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& getIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& getIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& getAdjointResidual() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& getAdjointResidual() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& getAdjointIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& getAdjointIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
