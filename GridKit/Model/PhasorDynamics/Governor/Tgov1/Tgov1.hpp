/**
 * @file Tgov1.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Wiktoria Zielinska (zielinskawa@ORNL.gov)
 * @brief Declaration of a Turbine Governor Model (IEEET1).
 *
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>

// Forward declarations
namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      template <typename RealT, typename IdxT>
      struct Tgov1Data;
    } // namespace Governor

    template <class ScalarT, typename IdxT>
    class Genrou;

    template <class ScalarT, typename IdxT>
    class SignalNode;

  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /// Internal variables of a `Tgov1`
      enum class Tgov1InternalVariables : size_t
      {
        PTX, ///< $P_{tx}$
        PV,  ///< $P_v$
        PM,  ///< $P_m$
        MAXIMUM,
      };

      /// External variables of a `Tgov1`
      enum class Tgov1ExternalVariables : size_t
      {
        DELTAOMEGA, ///< $\Delta_\omega$
        PREF,       ///< $P_{ref}$
        MAXIMUM,
      };

      template <class ScalarT, typename IdxT>
      class Tgov1 : public Component<ScalarT, IdxT>
      {
        using Component<ScalarT, IdxT>::gridkit_component_id_;
        using Component<ScalarT, IdxT>::alpha_;
        using Component<ScalarT, IdxT>::f_;
        using Component<ScalarT, IdxT>::nnz_;
        using Component<ScalarT, IdxT>::size_;
        using Component<ScalarT, IdxT>::tag_;
        using Component<ScalarT, IdxT>::time_;
        using Component<ScalarT, IdxT>::y_;
        using Component<ScalarT, IdxT>::yp_;
        using Component<ScalarT, IdxT>::wb_;
        using Component<ScalarT, IdxT>::h_;
        using Component<ScalarT, IdxT>::J_;

        using real_type       = typename Component<ScalarT, IdxT>::real_type;
        using model_data_type = Tgov1Data<real_type, IdxT>;
        using signal_type     = SignalNode<ScalarT, IdxT>;

      public:
        Tgov1();
        Tgov1(signal_type*, signal_type*);
        Tgov1(const model_data_type&);
        ~Tgov1() = default;

        int setGridKitComponentID(IdxT) override;
        int allocate() override;
        int verify() const override;
        int initialize() override;
        int tagDifferentiable() override;
        int evaluateResidual() override;

        // Still to be implemented
        int evaluateJacobian() override;

        void updateTime(real_type /* t */, real_type /* a */) override
        {
        }

        /// Get the `ComponentSignals` from this `Tgov1`
        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                Tgov1InternalVariables,
                                Tgov1ExternalVariables>&
        {
          return signals_;
        }

      public:
        __attribute__((always_inline)) inline int evaluateInternalResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);

      private:
        // Input parameters
        real_type R_{0};
        real_type Pvmin_{0};
        real_type Pvmax_{0};
        real_type T1_{0};
        real_type T2_{0};
        real_type T3_{0};
        real_type Dt_{0};

        // Input States (which can be parameters)
        ScalarT pref_{0};

        // Scale of Sigmoid function (temporary local implementation)
        const ScalarT mu_{4000.0};

        /// Component signal extension
        ComponentSignals<ScalarT, IdxT, Tgov1InternalVariables, Tgov1ExternalVariables> signals_;

        // Activation function (sigmoid approximation)
        ScalarT sigmoid(ScalarT x);

        // Indicator of Valve limit states
        ScalarT indicator_low(ScalarT x, ScalarT f);
        ScalarT indicator_high(ScalarT x, ScalarT f);
        ScalarT indicator(ScalarT x, ScalarT f);

        void initializeParameters(const model_data_type& data);

        /* Local copies of external variables */
        std::vector<ScalarT> ws_;
        std::map<IdxT, IdxT> ws_indices_;
      };

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
