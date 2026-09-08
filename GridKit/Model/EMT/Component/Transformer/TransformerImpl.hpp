#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>

#include <GridKit/Model/EMT/Component/Transformer/Transformer.hpp>
#include <GridKit/Model/EMT/Component/Transformer/TransformerData.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/EMT/PhasorInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Constructor for a three-phase transformer bank
     *
     * System sizes:
     * - Number of equations = 27
     * - Number of independent variables = 27
     */
    template <typename scalar_type, typename index_type>
    Transformer<scalar_type, index_type>::Transformer()
      : Transformer(ModelDataT{})
    {
    }

    template <typename scalar_type, typename index_type>
    Transformer<scalar_type, index_type>::Transformer(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      size_ = 27;
      setDerivedParams();
      for (size_t e = 0; e < 2; ++e)
      {
        for (size_t p = 0; p < 3; ++p)
        {
          assignOutput(static_cast<Outputs>(static_cast<size_t>(Outputs::i1a) + 3 * e + p), &current_[e][p]);
        }
      }
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    Transformer<scalar_type, index_type>::~Transformer()
    {
    }

    /**
     * @brief Read model parameters from the data object
     */
    template <typename scalar_type, typename index_type>
    void Transformer<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      S_              = parameter<RealT>(data, Parameter::S, S_);
      V1_             = parameter<RealT>(data, Parameter::V1, V1_);
      V2_             = parameter<RealT>(data, Parameter::V2, V2_);
      freq_           = parameter<RealT>(data, Parameter::f, freq_);
      P1_             = parameter<ABCMatrix<RealT>>(data, Parameter::P1, P1_);
      P2_             = parameter<ABCMatrix<RealT>>(data, Parameter::P2, P2_);
      tap_            = parameter<RealT>(data, Parameter::tap, tap_);
      R_              = parameter<RealT>(data, Parameter::R, R_);
      X_              = parameter<RealT>(data, Parameter::X, X_);
      I0_             = parameter<RealT>(data, Parameter::I0, I0_);
      P0_             = parameter<RealT>(data, Parameter::P0, P0_);
      knee_           = parameter<RealT>(data, Parameter::knee, knee_);
      Lsat_           = parameter<RealT>(data, Parameter::Lsat, TWO<RealT> * X_);
      split_          = parameter<RealT>(data, Parameter::split, split_);
    }

    /**
     * @brief Derived parameters
     *
     * The winding bases follow the positive-sequence winding voltage through
     * each connection map.
     */
    template <typename scalar_type, typename index_type>
    void Transformer<scalar_type, index_type>::setDerivedParams()
    {
      const RealT pi    = std::numbers::pi_v<RealT>;
      const RealT gamma = TWO<RealT> * pi / THREE<RealT>;

      omega_base_ = TWO<RealT> * pi * freq_;

      const std::array<std::complex<RealT>, 3> alpha{std::polar(ONE<RealT>, ZERO<RealT>),
                                                     std::polar(ONE<RealT>, -gamma),
                                                     std::polar(ONE<RealT>, gamma)};
      const std::array<RealT, 2>               rating{V1_, V2_};
      for (size_t e = 0; e < 2; ++e)
      {
        const auto&         map = connectionMap(e);
        std::complex<RealT> winding{};
        for (size_t n = 0; n < 3; ++n)
        {
          winding += map[n][0] * alpha[n];
        }
        const RealT v_winding = rating[e] / std::sqrt(THREE<RealT>) * std::abs(winding);
        v_peak_base_[e]       = std::numbers::sqrt2_v<RealT> * v_winding;
        i_peak_base_[e]       = 0.0;
        if (v_winding != 0.0)
        {
          i_peak_base_[e] = std::numbers::sqrt2_v<RealT> * S_ / (THREE<RealT> * v_winding);
        }
      }

      R1_ = HALF<RealT> * R_;
      R2_ = HALF<RealT> * R_;
      Gc_ = 0.0;
      if (S_ != 0.0)
      {
        Gc_ = P0_ / S_;
      }

      Lm_              = 0.0;
      inv_Lm_          = 0.0;
      const RealT i_m2 = I0_ * I0_ - Gc_ * Gc_;
      if (i_m2 > 0.0)
      {
        inv_Lm_ = std::sqrt(i_m2);
        Lm_     = ONE<RealT> / inv_Lm_;
      }
      k_sat_ = 0.0;
      if (Lsat_ != 0.0)
      {
        k_sat_ = ONE<RealT> / Lsat_ - inv_Lm_;
      }

      beta_[0] = split_;
      beta_[1] = ONE<RealT> - split_;
    }

    template <typename scalar_type, typename index_type>
    void Transformer<scalar_type, index_type>::attachTerminal(size_t end, PhaseSignals voltage)
    {
      if (allocated_)
        throw std::logic_error("Transformer terminals cannot change after allocation");
      if (end > 1)
        throw std::out_of_range("Invalid Transformer terminal");
      for (size_t p = 0; p < 3; ++p)
        if (!voltage[p])
          throw std::invalid_argument("Transformer requires terminal voltage inputs");
      for (size_t p = 0; p < 3; ++p)
        signals_.attachSignal(static_cast<TransformerExternalVariables>(3 * end + p), voltage[p]);
    }

    /**
     * @brief Bind a terminal current injection to an output signal.
     *
     * The injection folds the winding currents through the connection map
     * and scales them to SI amps.
     */
    template <typename scalar_type, typename index_type>
    void Transformer<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
    {
      if (output >= Outputs::SIZE || signal == nullptr)
        throw std::invalid_argument("Invalid Transformer output");
      const auto index = static_cast<size_t>(output) - static_cast<size_t>(Outputs::i1a);
      const auto end   = index / 3;
      const auto phase = index % 3;
      signal->claimProducer();
      signal->setComputed(
          [this, end, phase]
          { return terminalCurrent(end, phase); },
          [this, end, phase](typename SignalT::GradientT& gradient, RealT scale)
          {
            const auto& map = connectionMap(end);
            for (size_t k = 0; k < 3; ++k)
              gradient.emplace_back(this->getVariableIndex(static_cast<IdxT>(21 + 3 * end + k)), -scale * i_peak_base_[end] * map[phase][k]);
          });
    }

    template <typename scalar_type, typename index_type>
    scalar_type Transformer<scalar_type, index_type>::terminalCurrent(size_t end, size_t phase)
    {
      const auto* y     = y_.getData();
      const auto& map   = connectionMap(end);
      const auto  first = 21 + 3 * end;
      return -i_peak_base_[end] * (map[phase][0] * y[first] + map[phase][1] * y[first + 1] + map[phase][2] * y[first + 2]);
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method resizes local storage and registers coupling signals.
     */
    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }

      auto size = static_cast<size_t>(size_); // avoid compiler warnings

      tag_.resize(size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Resize coupling data
      this->allocateExternalVectors(static_cast<IdxT>(TransformerExternalVariables::MAXIMUM), 0);
      signals_.registerExternalVariableSignals(*this);
      signals_.bindInternalVariableSignals(*this);
      allocated_ = true;
      return 0;
    }

    /**
     * @brief Check model correctness
     *
     * @return Number of model configuration errors found
     */
    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::verify() const
    {
      int error_count = 0;

      if (!signals_.template isAttached<TransformerExternalVariables::V1A>()
          || !signals_.template isAttached<TransformerExternalVariables::V1B>()
          || !signals_.template isAttached<TransformerExternalVariables::V1C>()
          || !signals_.template isAttached<TransformerExternalVariables::V2A>()
          || !signals_.template isAttached<TransformerExternalVariables::V2B>()
          || !signals_.template isAttached<TransformerExternalVariables::V2C>())
      {
        Log::error() << "Transformer: a terminal voltage port is not attached\n";
        ++error_count;
      }

      if (S_ <= 0.0 || V1_ <= 0.0 || V2_ <= 0.0 || freq_ <= 0.0)
      {
        Log::error() << "Transformer: the ratings must be positive\n";
        ++error_count;
      }

      if (tap_ <= 0.0)
      {
        Log::error() << "Transformer: the tap ratio must be positive\n";
        ++error_count;
      }

      if (R_ < 0.0 || P0_ < 0.0)
      {
        Log::error() << "Transformer: the short-circuit resistance and no-load "
                        "loss must be nonnegative\n";
        ++error_count;
      }

      if (X_ <= 0.0 || knee_ <= 0.0 || Lsat_ <= 0.0)
      {
        Log::error() << "Transformer: the short-circuit reactance, knee flux "
                        "linkage, and saturation inductance must be positive\n";
        ++error_count;
      }

      if (I0_ * S_ <= P0_)
      {
        Log::error() << "Transformer: the no-load current must exceed the "
                        "no-load loss share\n";
        ++error_count;
      }

      if (split_ < 0.0 || split_ > 1.0)
      {
        Log::error() << "Transformer: the magnetizing share must lie in [0, 1]\n";
        ++error_count;
      }

      for (size_t e = 0; e < 2; ++e)
      {
        const auto& map       = connectionMap(e);
        bool        map_valid = true;
        for (size_t k = 0; k < 3; ++k)
        {
          int plus  = 0;
          int minus = 0;
          for (size_t n = 0; n < 3; ++n)
          {
            if (map[n][k] == 1.0)
            {
              ++plus;
            }
            else if (map[n][k] == -1.0)
            {
              ++minus;
            }
            else if (map[n][k] != 0.0)
            {
              map_valid = false;
            }
          }
          if (plus != 1 || minus > 1)
          {
            map_valid = false;
          }
        }
        if (!map_valid)
        {
          Log::error() << "Transformer: each connection map column needs one +1, "
                          "at most one -1, and zeros elsewhere\n";
          ++error_count;
        }
      }

      if (Lm_ <= Lsat_)
      {
        Log::error() << "Transformer: the magnetizing inductance must exceed the "
                        "saturation inductance\n";
        ++error_count;
      }

      return error_count;
    }

    /**
     * Initialization of the transformer model
     *
     * The series current and flux linkages start de-energized unless flux
     * linkages are supplied. The assembled consistent initialization resolves
     * the algebraic variables.
     */
    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::initialize(const std::array<RealT, 6>& flux)
    {
      for (const auto value : flux)
      {
        if (!std::isfinite(value))
          throw std::invalid_argument("Transformer: nonfinite flux linkage state");
      }
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      for (IdxT j = 0; j < size_; ++j)
      {
        y[j]  = 0.0;
        yp[j] = 0.0;
      }
      for (size_t n = 0; n < 6; ++n)
      {
        y[3 + n] = static_cast<ScalarT>(flux[n]);
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Transformer<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
    {
      for (const auto& [key, value] : values)
        if (std::find(flux_keys_.begin(), flux_keys_.end(), key) == flux_keys_.end() || !std::isfinite(value))
          throw std::invalid_argument("Transformer: invalid initial state " + key);
    }

    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
    {
      validateInitialState(values);
      std::array<RealT, 6> flux{};
      for (size_t n = 0; n < 6; ++n)
      {
        const auto entry = values.find(std::string(flux_keys_[n]));
        if (entry != values.end())
          flux[n] = entry->second;
      }
      return initialize(flux);
    }

    /**
     * @brief Solve the unsaturated phasor circuit from the terminal samples.
     *
     * The node voltages of each phase satisfy a two-by-two complex system,
     * embedded in the three-phase phasor solver with an identity third row.
     * The remaining states follow from the sinusoidal integrals of the node
     * voltages.
     */
    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::initializeSteadyState(RealT omega)
    {
      if (!std::isfinite(omega) || omega <= ZERO<RealT>)
      {
        return 1;
      }
      if (omega_base_ <= ZERO<RealT> || X_ <= ZERO<RealT> || Lm_ <= ZERO<RealT>
          || v_peak_base_[0] <= ZERO<RealT> || v_peak_base_[1] <= ZERO<RealT>)
      {
        return 1;
      }
      this->gatherExternalVariables();

      const RealT      w   = omega / omega_base_;
      const RealT      wx  = w * X_;
      const RealT      wlm = w * Lm_;
      const RealT      t2  = tap_ * tap_;
      ABCMatrix<RealT> H_re{};
      ABCMatrix<RealT> H_im{};
      H_re[0][0] = ONE<RealT> + R1_ * beta_[0] * Gc_;
      H_im[0][0] = -R1_ * beta_[0] / wlm - R1_ / wx;
      H_im[0][1] = R1_ / wx;
      H_im[1][0] = t2 * R2_ / wx;
      H_re[1][1] = ONE<RealT> + t2 * R2_ * beta_[1] * Gc_;
      H_im[1][1] = -t2 * R2_ * beta_[1] / wlm - t2 * R2_ / wx;
      H_re[2][2] = ONE<RealT>;

      auto* y  = y_.getData();
      auto* yp = yp_.getData();
      for (size_t n = 0; n < 3; ++n)
      {
        // Winding voltage samples in per unit through the connection maps
        ABCVector<RealT> v{};
        ABCVector<RealT> v_dot{};
        for (size_t m = 0; m < 3; ++m)
        {
          v[0]     += P1_[m][n] * static_cast<RealT>(y_ext_[m]) / v_peak_base_[0];
          v_dot[0] += P1_[m][n] * static_cast<RealT>(yp_ext_[m]) / v_peak_base_[0];
          v[1]     += tap_ * P2_[m][n] * static_cast<RealT>(y_ext_[3 + m]) / v_peak_base_[1];
          v_dot[1] += tap_ * P2_[m][n] * static_cast<RealT>(yp_ext_[3 + m]) / v_peak_base_[1];
        }
        ABCVector<RealT> e{};
        ABCVector<RealT> e_dot{};
        if (solvePhasorSystem(H_re, H_im, omega, v, v_dot, e, e_dot) != 0)
        {
          return 1;
        }

        // Sinusoidal integrals of the node voltages
        const RealT psi1     = -omega_base_ * e_dot[0] / (omega * omega);
        const RealT psi2     = -omega_base_ * e_dot[1] / (omega * omega);
        const RealT psi1_dot = omega_base_ * e[0];
        const RealT psi2_dot = omega_base_ * e[1];
        const RealT i12      = -omega_base_ * (e_dot[0] - e_dot[1]) / (X_ * omega * omega);
        const RealT i12_dot  = omega_base_ * (e[0] - e[1]) / X_;
        const RealT im1      = beta_[0] * (inv_Lm_ * psi1 + Gc_ * e[0]);
        const RealT im2      = beta_[1] * (inv_Lm_ * psi2 + Gc_ * e[1]);
        const RealT im1_dot  = beta_[0] * (inv_Lm_ * psi1_dot + Gc_ * e_dot[0]);
        const RealT im2_dot  = beta_[1] * (inv_Lm_ * psi2_dot + Gc_ * e_dot[1]);

        y[n]       = i12;
        y[3 + n]   = psi1;
        y[6 + n]   = psi2;
        y[9 + n]   = e[0];
        y[12 + n]  = e[1];
        y[15 + n]  = im1;
        y[18 + n]  = im2;
        y[21 + n]  = im1 + i12;
        y[24 + n]  = tap_ * (im2 - i12);
        yp[n]      = i12_dot;
        yp[3 + n]  = psi1_dot;
        yp[6 + n]  = psi2_dot;
        yp[9 + n]  = e_dot[0];
        yp[12 + n] = e_dot[1];
        yp[15 + n] = im1_dot;
        yp[18 + n] = im2_dot;
        yp[21 + n] = im1_dot + i12_dot;
        yp[24 + n] = tap_ * (im2_dot - i12_dot);
      }
      y_.setDataUpdated();
      yp_.setDataUpdated();
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam scalar_type Scalar data type
     * @tparam index_type Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline scalar_type Transformer<scalar_type, index_type>::magnetizingCurrent(ScalarT psi) const
    {
      return inv_Lm_ * psi + k_sat_ * (Math::ramp(psi - knee_) - Math::ramp(-psi - knee_));
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Transformer<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        const ScalarT*                  yp,
        const ScalarT*                  y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f)
    {
      /* Read variables */
      const ScalarT i12a  = y[0];
      const ScalarT i12b  = y[1];
      const ScalarT i12c  = y[2];
      const ScalarT psi1a = y[3];
      const ScalarT psi1b = y[4];
      const ScalarT psi1c = y[5];
      const ScalarT psi2a = y[6];
      const ScalarT psi2b = y[7];
      const ScalarT psi2c = y[8];
      const ScalarT e1a   = y[9];
      const ScalarT e1b   = y[10];
      const ScalarT e1c   = y[11];
      const ScalarT e2a   = y[12];
      const ScalarT e2b   = y[13];
      const ScalarT e2c   = y[14];
      const ScalarT im1a  = y[15];
      const ScalarT im1b  = y[16];
      const ScalarT im1c  = y[17];
      const ScalarT im2a  = y[18];
      const ScalarT im2b  = y[19];
      const ScalarT im2c  = y[20];
      const ScalarT iw1a  = y[21];
      const ScalarT iw1b  = y[22];
      const ScalarT iw1c  = y[23];
      const ScalarT iw2a  = y[24];
      const ScalarT iw2b  = y[25];
      const ScalarT iw2c  = y[26];

      /* Read derivatives */
      const ScalarT i12a_dot  = yp[0];
      const ScalarT i12b_dot  = yp[1];
      const ScalarT i12c_dot  = yp[2];
      const ScalarT psi1a_dot = yp[3];
      const ScalarT psi1b_dot = yp[4];
      const ScalarT psi1c_dot = yp[5];
      const ScalarT psi2a_dot = yp[6];
      const ScalarT psi2b_dot = yp[7];
      const ScalarT psi2c_dot = yp[8];

      // Set coupling variable aliases
      const ScalarT v1a = y_ext[0];
      const ScalarT v1b = y_ext[1];
      const ScalarT v1c = y_ext[2];
      const ScalarT v2a = y_ext[3];
      const ScalarT v2b = y_ext[4];
      const ScalarT v2c = y_ext[5];

      // Winding voltages in transformer per unit through the connection maps
      const ScalarT vw1a = (P1_[0][0] * v1a + P1_[1][0] * v1b + P1_[2][0] * v1c) / v_peak_base_[0];
      const ScalarT vw1b = (P1_[0][1] * v1a + P1_[1][1] * v1b + P1_[2][1] * v1c) / v_peak_base_[0];
      const ScalarT vw1c = (P1_[0][2] * v1a + P1_[1][2] * v1b + P1_[2][2] * v1c) / v_peak_base_[0];
      const ScalarT vw2a = (P2_[0][0] * v2a + P2_[1][0] * v2b + P2_[2][0] * v2c) / v_peak_base_[1];
      const ScalarT vw2b = (P2_[0][1] * v2a + P2_[1][1] * v2b + P2_[2][1] * v2c) / v_peak_base_[1];
      const ScalarT vw2c = (P2_[0][2] * v2a + P2_[1][2] * v2b + P2_[2][2] * v2c) / v_peak_base_[1];

      const RealT inv_omega_base = ONE<RealT> / omega_base_;
      const RealT x_over_omega   = X_ * inv_omega_base;

      /* 9 transformer differential equations */
      f[0] = x_over_omega * i12a_dot + e2a - e1a;
      f[1] = x_over_omega * i12b_dot + e2b - e1b;
      f[2] = x_over_omega * i12c_dot + e2c - e1c;
      f[3] = inv_omega_base * psi1a_dot - e1a;
      f[4] = inv_omega_base * psi1b_dot - e1b;
      f[5] = inv_omega_base * psi1c_dot - e1c;
      f[6] = inv_omega_base * psi2a_dot - e2a;
      f[7] = inv_omega_base * psi2b_dot - e2b;
      f[8] = inv_omega_base * psi2c_dot - e2c;

      /* 18 transformer algebraic equations */
      f[9]  = e1a - vw1a + R1_ * iw1a;
      f[10] = e1b - vw1b + R1_ * iw1b;
      f[11] = e1c - vw1c + R1_ * iw1c;
      f[12] = e2a - tap_ * vw2a + tap_ * R2_ * iw2a;
      f[13] = e2b - tap_ * vw2b + tap_ * R2_ * iw2b;
      f[14] = e2c - tap_ * vw2c + tap_ * R2_ * iw2c;
      f[15] = im1a - beta_[0] * magnetizingCurrent(psi1a) - beta_[0] * Gc_ * e1a;
      f[16] = im1b - beta_[0] * magnetizingCurrent(psi1b) - beta_[0] * Gc_ * e1b;
      f[17] = im1c - beta_[0] * magnetizingCurrent(psi1c) - beta_[0] * Gc_ * e1c;
      f[18] = im2a - beta_[1] * magnetizingCurrent(psi2a) - beta_[1] * Gc_ * e2a;
      f[19] = im2b - beta_[1] * magnetizingCurrent(psi2b) - beta_[1] * Gc_ * e2b;
      f[20] = im2c - beta_[1] * magnetizingCurrent(psi2c) - beta_[1] * Gc_ * e2c;
      f[21] = iw1a - im1a - i12a;
      f[22] = iw1b - im1b - i12b;
      f[23] = iw1c - im1c - i12c;
      f[24] = iw2a + tap_ * i12a - tap_ * im2a;
      f[25] = iw2b + tap_ * i12b - tap_ * im2b;
      f[26] = iw2c + tap_ * i12c - tap_ * im2c;

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::evaluateInternalResidual()
    {
      this->gatherExternalVariables();

      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      auto*       f  = f_.getData();
      evaluateInternalResidual(y, yp, y_ext_.data(), yp_ext_.data(), f);
      f_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Assemble the transformer equations.
     *
     */
    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::evaluateResidual()
    {
      return evaluateInternalResidual();
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Transformer<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Transformer<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;

      monitor_->set(Variable::i12a, [this]
                    { return y_.getData()[0]; });
      monitor_->set(Variable::i12b, [this]
                    { return y_.getData()[1]; });
      monitor_->set(Variable::i12c, [this]
                    { return y_.getData()[2]; });
      monitor_->set(Variable::psi1a, [this]
                    { return y_.getData()[3]; });
      monitor_->set(Variable::psi1b, [this]
                    { return y_.getData()[4]; });
      monitor_->set(Variable::psi1c, [this]
                    { return y_.getData()[5]; });
      monitor_->set(Variable::psi2a, [this]
                    { return y_.getData()[6]; });
      monitor_->set(Variable::psi2b, [this]
                    { return y_.getData()[7]; });
      monitor_->set(Variable::psi2c, [this]
                    { return y_.getData()[8]; });
      monitor_->set(Variable::i1a, [this]
                    { return terminalCurrent(0, 0); });
      monitor_->set(Variable::i1b, [this]
                    { return terminalCurrent(0, 1); });
      monitor_->set(Variable::i1c, [this]
                    { return terminalCurrent(0, 2); });
      monitor_->set(Variable::i2a, [this]
                    { return terminalCurrent(1, 0); });
      monitor_->set(Variable::i2b, [this]
                    { return terminalCurrent(1, 1); });
      monitor_->set(Variable::i2c, [this]
                    { return terminalCurrent(1, 2); });
    }

  } // namespace EMT
} // namespace GridKit
