# PhasorDynamics `Container`: two-stage implementation plan

## Summary

Implement `Container<ScalarT, IdxT> : Component<ScalarT, IdxT>` as a component-only composite:

- A concrete composite owns its child components as members; `Container` keeps non-owning pointers.
- The same `Container` can be nested under another `Container`, added to `SystemModel`, or used directly as a solver root.
- Children share slices of one root `y_`, `yp_`, `f_`, and `abs_tol_`; no copying and no `Vector` changes.
- Public input/output ports are ordinary `SignalNode`s forwarded unchanged to children.
- “Internal” means owned by the component tree; “external” means owned elsewhere. It does not mean private versus public.
- Internal residuals are assigned globally before any external residual contributions are added.
- A direct root must be closed: all Jacobian rows and columns must belong to that root.
- Topology is assembled before allocation and fixed afterward.
- No buses, monitors, JSON/factory support, new port vocabulary, or migration of existing model outputs are included.

Use two stacked review units:

1. Connection and residual-phase foundation.
2. Header-only `Container` plus its dedicated test target.

The implementation is limited to the code shapes below. Surrounding existing code remains unchanged.

## 1. Target usage

A composite exposes its ports through the current `ComponentSignals` API. Its output can be public while remaining internally owned:

```cpp
enum class CascadeInternalVariables : size_t
{
  X,
  MAXIMUM
};

enum class CascadeExternalVariables : size_t
{
  U,
  MAXIMUM
};

template <typename ScalarT, typename IdxT>
class Cascade : public PhasorDynamics::Container<ScalarT, IdxT>
{
public:
  using ContainerT = PhasorDynamics::Container<ScalarT, IdxT>;
  using SignalsT   = PhasorDynamics::ComponentSignals<
      ScalarT,
      IdxT,
      CascadeInternalVariables,
      CascadeExternalVariables>;

  Cascade()
    : first_(3.0, 4.0, 5.0, 3.0, 0.2),
      second_(6.0, 7.0, 8.0, 5.0, 0.3)
  {
    first_.getSignals()
        .template assignSignalNode<LinearInternalVariables::X>(&middle_);
    second_.getSignals()
        .template attachSignalNode<LinearExternalVariables::U>(&middle_);

    this->addComponent(&first_);
    this->addComponent(&second_);
  }

  int allocate() override
  {
    if (signals_.template isAttached<CascadeExternalVariables::U>())
    {
      first_.getSignals()
          .template attachSignalNode<LinearExternalVariables::U>(
              signals_.template getSignalNode<CascadeExternalVariables::U>());
    }

    if (signals_.template isAssigned<CascadeInternalVariables::X>())
    {
      second_.getSignals()
          .template assignSignalNode<LinearInternalVariables::X>(
              signals_.template getSignalNode<CascadeInternalVariables::X>());
    }

    return ContainerT::allocate();
  }

  auto getSignals() -> SignalsT&
  {
    return signals_;
  }

private:
  SignalsT                          signals_;
  PhasorDynamics::SignalNode<ScalarT, IdxT> middle_;
  LinearComponent<ScalarT, IdxT>   first_;
  LinearComponent<ScalarT, IdxT>   second_;
};
```

A closed direct-solver root wires the composite to another component with the exact same nodes:

```cpp
template <typename ScalarT, typename IdxT>
class ClosedModel : public PhasorDynamics::Container<ScalarT, IdxT>
{
public:
  ClosedModel()
    : boundary_(1.0, 9.0, 2.0, 2.0, 0.1)
  {
    boundary_.getSignals()
        .template assignSignalNode<LinearInternalVariables::X>(&u_);
    boundary_.getSignals()
        .template attachSignalNode<LinearExternalVariables::U>(&z_);

    cascade_.getSignals()
        .template attachSignalNode<CascadeExternalVariables::U>(&u_);
    cascade_.getSignals()
        .template assignSignalNode<CascadeInternalVariables::X>(&z_);

    this->addComponent(&boundary_);
    this->addComponent(&cascade_);
  }

  auto uNode() -> PhasorDynamics::SignalNode<ScalarT, IdxT>*
  {
    return &u_;
  }

  auto zNode() -> PhasorDynamics::SignalNode<ScalarT, IdxT>*
  {
    return &z_;
  }

  auto cascade() -> Cascade<ScalarT, IdxT>&
  {
    return cascade_;
  }

private:
  PhasorDynamics::SignalNode<ScalarT, IdxT> u_;
  PhasorDynamics::SignalNode<ScalarT, IdxT> z_;
  LinearComponent<ScalarT, IdxT>            boundary_;
  Cascade<ScalarT, IdxT>                    cascade_;
};
```

It is directly usable through the existing evaluator contract:

```cpp
ClosedModel<double, size_t> model;

GridKit::Model::Evaluator<double, size_t>* evaluator = &model;

evaluator->allocate();
evaluator->initialize();
evaluator->updateTime(0.0, 10.0);
evaluator->tagDifferentiable();
evaluator->setAbsoluteTolerance(1.0e-6);
evaluator->evaluateResidual();
evaluator->evaluateJacobian();
```

The storage contract is:

| Ownership | Values | Derivatives | Residuals | Indices |
|---|---|---|---|---|
| Owned/internal | `y_int_` | `yp_int_` | `f_int_` | existing `variable_indices_`, `residual_indices_` |
| External aliases | `y_ext_` | `yp_ext_` | `f_ext_` | `variable_indices_ext_`, `residual_indices_ext_` |

External residual writes are always additive:

```cpp
*f_ext_[i] += contribution;
```

No `abs_tol_ext_` or `tag_ext_` is needed because tolerance and differential tags belong to the owner.

## 2. Stage one: connection foundation

### `SignalNode`

Keep the current two-pointer overload for constants and adapters. Add the complete solver-variable overload without changing `linked()` semantics.

```cpp
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class Component;

    template <typename scalar_type, typename index_type>
    class SignalNode
    {
      template <typename, typename>
      friend class Component;

    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      SignalNode();
      SignalNode(const SignalNodeData<RealT, IdxT>& data);

      virtual ~SignalNode() = default;

      void set(ScalarT* signal, IdxT* variable_index);

      void set(ScalarT* signal,
               ScalarT* yp,
               ScalarT* f,
               IdxT*    variable_index,
               IdxT*    residual_index);

      bool    linked() const;
      ScalarT read() const;
      void    init(ScalarT signal);

      const IdxT signalId() const
      {
        return signal_id_;
      }

      IdxT getVariableIndex() const
      {
        return *variable_index_;
      }

      virtual const IdxT busID() const
      {
        return bus_id_;
      }

    private:
      ScalarT* signal_{nullptr};
      ScalarT* yp_{nullptr};
      ScalarT* f_{nullptr};
      IdxT*    residual_index_{nullptr};
      IdxT     signal_id_{0};

    protected:
      const IdxT bus_id_{INVALID_INDEX<IdxT>};
      IdxT*      variable_index_{nullptr};
    };
  }
}
```

Implementation:

```cpp
template <typename scalar_type, typename index_type>
void SignalNode<scalar_type, index_type>::set(
    ScalarT* signal,
    IdxT*    variable_index)
{
  set(signal, nullptr, nullptr, variable_index, nullptr);
}

template <typename scalar_type, typename index_type>
void SignalNode<scalar_type, index_type>::set(
    ScalarT* signal,
    ScalarT* yp,
    ScalarT* f,
    IdxT*    variable_index,
    IdxT*    residual_index)
{
  signal_         = signal;
  yp_             = yp;
  f_              = f;
  variable_index_ = variable_index;
  residual_index_ = residual_index;
}

template <typename scalar_type, typename index_type>
bool SignalNode<scalar_type, index_type>::linked() const
{
  return signal_ != nullptr && variable_index_ != nullptr;
}
```

Delegating the legacy overload clears any previous `yp_`, `f_`, and residual-index aliases. Existing value-only signals remain valid.

### `Component`

Add the existing `SignalNode` type and make only `size()` overridable; `nnz()` remains final.

```cpp
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>

// ...

using SignalT = SignalNode<ScalarT, IdxT>;

IdxT size() override
{
  return size_;
}
```

Replace the vector-slice `bind` body with a zero-state-safe version. `bound_` is set only by a successful vector-slice binding:

```cpp
int bind(
    VectorT& y,
    VectorT& yp,
    VectorT& f,
    VectorT& abs_tol,
    IdxT    offset)
{
  const IdxT n = size();

  if (offset > y.getSize() || n > y.getSize() - offset
      || offset > yp.getSize() || n > yp.getSize() - offset
      || offset > f.getSize() || n > f.getSize() - offset
      || offset > abs_tol.getSize() || n > abs_tol.getSize() - offset)
  {
    Log::error() << "Component::bind - system vectors are smaller than "
                 << "offset + size = " << offset + n << "\n";
    return 1;
  }

  ScalarT* y_data       = nullptr;
  ScalarT* yp_data      = nullptr;
  ScalarT* f_data       = nullptr;
  ScalarT* abs_tol_data = nullptr;

  if (n != IdxT{})
  {
    auto* system_y       = y.getData(memory::HOST);
    auto* system_yp      = yp.getData(memory::HOST);
    auto* system_f       = f.getData(memory::HOST);
    auto* system_abs_tol = abs_tol.getData(memory::HOST);

    if (system_y == nullptr || system_yp == nullptr
        || system_f == nullptr || system_abs_tol == nullptr)
    {
      Log::error() << "Component::bind - system vector data is null or stale\n";
      return 1;
    }

    y_data       = system_y + offset;
    yp_data      = system_yp + offset;
    f_data       = system_f + offset;
    abs_tol_data = system_abs_tol + offset;
  }

  const int y_status =
      y_.setData(y_data, n, memory::HOST);
  const int yp_status =
      yp_.setData(yp_data, n, memory::HOST);
  const int f_status =
      f_.setData(f_data, n, memory::HOST);
  const int abs_tol_status =
      abs_tol_.setData(abs_tol_data, n, memory::HOST);

  if (y_status != 0 || yp_status != 0
      || f_status != 0 || abs_tol_status != 0)
  {
    Log::error() << "Component::bind - failed to bind vectors to system storage\n";
    return 1;
  }

  y_int_   = y_data;
  yp_int_  = yp_data;
  f_int_   = f_data;
  bound_   = true;
  allocated_ = true;

  return 0;
}
```

Add the overloaded signal binding. It records topology but deliberately does not resolve the pointers yet:

```cpp
int bind(SignalT* signal, IdxT external_index)
{
  if (signal == nullptr || external_index == INVALID_INDEX<IdxT>)
  {
    Log::error() << "Component::bind - invalid external signal\n";
    return 1;
  }

  const auto slot = static_cast<size_t>(external_index);
  const auto count = slot + 1;

  signal_nodes_ext_.resize(count, nullptr);
  y_ext_.resize(count, nullptr);
  yp_ext_.resize(count, nullptr);
  f_ext_.resize(count, nullptr);
  variable_indices_ext_.resize(count, nullptr);
  residual_indices_ext_.resize(count, nullptr);

  signal_nodes_ext_[slot]       = signal;
  y_ext_[slot]                  = nullptr;
  yp_ext_[slot]                 = nullptr;
  f_ext_[slot]                  = nullptr;
  variable_indices_ext_[slot]   = nullptr;
  residual_indices_ext_[slot]   = nullptr;

  return 0;
}
```

Resolve recorded nodes only after all producers have allocated and global indices have been assigned:

```cpp
virtual int bindSignalNodes(IdxT system_size)
{
  for (size_t i = 0; i < signal_nodes_ext_.size(); ++i)
  {
    auto* signal = signal_nodes_ext_[i];

    if (signal == nullptr)
    {
      continue;
    }

    if (signal->signal_ == nullptr
        || signal->yp_ == nullptr
        || signal->f_ == nullptr
        || signal->variable_index_ == nullptr
        || signal->residual_index_ == nullptr)
    {
      Log::error() << "Component::bindSignalNodes - external signal "
                   << i << " is not a complete solver-variable signal\n";
      return 1;
    }

    const IdxT variable_index = *signal->variable_index_;
    const IdxT residual_index = *signal->residual_index_;

    if (variable_index == INVALID_INDEX<IdxT>
        || residual_index == INVALID_INDEX<IdxT>
        || variable_index >= system_size
        || residual_index >= system_size)
    {
      Log::error() << "Component::bindSignalNodes - external signal "
                   << i << " is outside the root system\n";
      return 1;
    }

    y_ext_[i]                = signal->signal_;
    yp_ext_[i]               = signal->yp_;
    f_ext_[i]                = signal->f_;
    variable_indices_ext_[i] = signal->variable_index_;
    residual_indices_ext_[i] = signal->residual_index_;
  }

  return 0;
}
```

Add compatibility residual phases. Existing components require no edits:

```cpp
virtual int evaluateInternalResidual()
{
  return evaluateResidual();
}

virtual int evaluateExternalResidual()
{
  return 0;
}
```

Make the existing propagation methods virtual and guard their indices:

```cpp
virtual int setVariableIndex(IdxT local_index, IdxT global_index)
{
  if (local_index >= static_cast<IdxT>(variable_indices_.size()))
  {
    Log::error() << "Component::setVariableIndex - local index is out of range\n";
    return 1;
  }

  variable_indices_[static_cast<size_t>(local_index)] = global_index;
  return 0;
}

virtual int setResidualIndex(IdxT local_index, IdxT global_index)
{
  if (local_index >= static_cast<IdxT>(residual_indices_.size()))
  {
    Log::error() << "Component::setResidualIndex - local index is out of range\n";
    return 1;
  }

  residual_indices_[static_cast<size_t>(local_index)] = global_index;
  return 0;
}

virtual void setSystemBase(
    RealT freq_system_base,
    RealT va_system_base)
{
  freq_system_base_ = freq_system_base;
  va_system_base_   = va_system_base;
}
```

Refresh the internal aliases whenever local storage is allocated:

```cpp
void allocateVectors(IdxT n)
{
  y_.resize(n);
  yp_.resize(n);
  f_.resize(n);
  abs_tol_.resize(n);

  y_int_  = y_.getData(memory::HOST);
  yp_int_ = yp_.getData(memory::HOST);
  f_int_  = f_.getData(memory::HOST);
}
```

Add the protected storage:

```cpp
ScalarT* y_int_{nullptr};
ScalarT* yp_int_{nullptr};
ScalarT* f_int_{nullptr};

std::vector<const ScalarT*> y_ext_;
std::vector<const ScalarT*> yp_ext_;
std::vector<ScalarT*>       f_ext_;

std::vector<const IdxT*> variable_indices_ext_;
std::vector<const IdxT*> residual_indices_ext_;
std::vector<SignalT*>    signal_nodes_ext_;

bool bound_{false};
```

### `ComponentSignals`

Add an external-variable overload of the existing name. The constraint prevents redeclaration when both enum types are identical:

```cpp
template <ExternalVariables variable>
  requires(!std::is_same_v<InternalVariables, ExternalVariables>)
auto getSignalNode() -> SignalNode<ScalarT, IdxT>*
{
  static_assert(variable < ExternalVariables::MAXIMUM);

  if (!external_variable_signals_[static_cast<size_t>(variable)])
  {
    throw std::logic_error(
        "A signal node has not been attached to this external variable");
  }

  return *external_variable_signals_[static_cast<size_t>(variable)];
}
```

### `SystemModel`

Add these declarations:

```cpp
int bindSignalNodes(IdxT system_size) override;

int evaluateInternalResidual() override;
int evaluateExternalResidual() override;

void setSystemBase(
    RealT freq_system_base,
    RealT va_system_base) override;
```

After all bus/component allocation and global-index loops, but before `verify()`, add:

```cpp
if (this->bindSignalNodes(size_) != 0)
{
  Log::error() << "Failed to bind component signal nodes\n";
  throw std::runtime_error("SystemModel allocation failed");
}
```

Before the existing allocation-time sparse evaluation, initialize time and `alpha_`:

```cpp
if (hasJacobian())
{
  updateTime(RealT{0}, RealT{1});
  initialize();
  evaluateResidual();
  evaluateJacobian();
}
```

Implement recursive post-allocation binding:

```cpp
template <typename scalar_type, typename index_type>
int SystemModel<scalar_type, index_type>::bindSignalNodes(IdxT system_size)
{
  int status = ComponentT::bindSignalNodes(system_size);

  for (const auto& component : components_)
  {
    if (component->bindSignalNodes(system_size) != 0)
    {
      status = 1;
    }
  }

  return status;
}
```

Replace the residual traversal with the two global phases:

```cpp
template <typename scalar_type, typename index_type>
int SystemModel<scalar_type, index_type>::evaluateResidual()
{
  evaluateInternalResidual();
  evaluateExternalResidual();

  f_.setDataUpdated();
  return 0;
}

template <typename scalar_type, typename index_type>
int SystemModel<scalar_type, index_type>::evaluateInternalResidual()
{
  for (const auto& bus : buses_)
  {
    bus->evaluateResidual();
  }

  for (const auto& component : components_)
  {
    component->evaluateInternalResidual();
  }

  return 0;
}

template <typename scalar_type, typename index_type>
int SystemModel<scalar_type, index_type>::evaluateExternalResidual()
{
  for (const auto& component : components_)
  {
    component->evaluateExternalResidual();
  }

  return 0;
}
```

Buses remain first because they assign bus residuals before components add to them. Child return codes continue to be ignored, matching current `SystemModel` traversal behavior.

### Stage-one focused tests

Add a zero-state phase probe to `SystemTests.hpp`:

```cpp
template <typename ScalarT, typename IdxT>
class ResidualPhaseProbe
  : public PhasorDynamics::Component<ScalarT, IdxT>
{
public:
  using ComponentT = PhasorDynamics::Component<ScalarT, IdxT>;
  using RealT      = typename ComponentT::RealT;

  ResidualPhaseProbe(std::vector<int>& calls, int id)
    : calls_(calls),
      id_(id)
  {
  }

  int setGridKitComponentID(IdxT id) override
  {
    this->gridkit_component_id_ = id;
    return 0;
  }

  int allocate() override
  {
    this->allocated_ = true;
    return 0;
  }

  int verify() const override
  {
    return 0;
  }

  int initialize() override
  {
    return 0;
  }

  int tagDifferentiable() override
  {
    return 0;
  }

  int setAbsoluteTolerance(RealT) override
  {
    return 0;
  }

  int evaluateResidual() override
  {
    evaluateInternalResidual();
    evaluateExternalResidual();
    return 0;
  }

  int evaluateInternalResidual() override
  {
    calls_.push_back(id_);
    return 0;
  }

  int evaluateExternalResidual() override
  {
    calls_.push_back(-id_);
    return 0;
  }

  int evaluateJacobian() override
  {
    return 0;
  }

  bool hasJacobian() override
  {
    return false;
  }

private:
  std::vector<int>& calls_;
  int               id_;
};
```

Add these test methods:

```cpp
TestOutcome signalNodeOverloads()
{
  TestStatus success = true;

  ScalarT y{2.0};
  ScalarT yp{3.0};
  ScalarT f{4.0};
  IdxT    variable_index{5};
  IdxT    residual_index{6};

  PhasorDynamics::SignalNode<ScalarT, IdxT> signal;

  signal.set(&y, &variable_index);
  success *= signal.linked();
  success *= isEqual(signal.read(), ScalarT{2.0});

  signal.init(ScalarT{7.0});
  success *= isEqual(y, ScalarT{7.0});

  signal.set(
      &y,
      &yp,
      &f,
      &variable_index,
      &residual_index);

  success *= signal.linked();
  success *= signal.getVariableIndex() == variable_index;

  return success.report(__func__);
}

TestOutcome residualPhases()
{
  TestStatus success = true;

  std::vector<int> calls;

  ResidualPhaseProbe<ScalarT, IdxT> first(calls, 1);
  ResidualPhaseProbe<ScalarT, IdxT> second(calls, 2);

  PhasorDynamics::SystemModel<ScalarT, IdxT> system;
  system.addComponent(&first);
  system.addComponent(&second);

  success *= system.allocate() == 0;
  success *= system.evaluateResidual() == 0;
  success *= calls == std::vector<int>({1, 2, -1, -2});

  return success.report(__func__);
}
```

Add the corresponding runner calls:

```cpp
result += test.signalNodeOverloads();
result += test.residualPhases();
```

No stage-one CMake changes are required.

## 3. Stage two: `Container`

### `Container.hpp`

```cpp
#pragma once

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class Container : public Component<scalar_type, index_type>
    {
    public:
      using ComponentT = Component<scalar_type, index_type>;
      using ScalarT    = typename ComponentT::ScalarT;
      using IdxT       = typename ComponentT::IdxT;
      using RealT      = typename ComponentT::RealT;
      using CooMatrixT = typename ComponentT::CooMatrixT;

      Container() = default;
      ~Container() override = default;

      Container(const Container&) = delete;
      Container& operator=(const Container&) = delete;
      Container(Container&&) = delete;
      Container& operator=(Container&&) = delete;

      IdxT size() override;

      int setGridKitComponentID(IdxT component_id) override;

      int allocate() override;
      int verify() const override;
      int initialize() override;
      int bindSignalNodes(IdxT system_size) override;

      int tagDifferentiable() override;
      int setAbsoluteTolerance(RealT rel_tol) override;

      int evaluateResidual() override;
      int evaluateInternalResidual() override;
      int evaluateExternalResidual() override;
      int evaluateJacobian() override;

      bool hasJacobian() override;

      void updateTime(RealT t, RealT alpha) override;
      void setSystemBase(
          RealT freq_system_base,
          RealT va_system_base) override;

      int setVariableIndex(
          IdxT local_index,
          IdxT global_index) override;

      int setResidualIndex(
          IdxT local_index,
          IdxT global_index) override;

      void addComponent(ComponentT* component);
      ComponentT* getComponent(IdxT component_id);

    private:
      std::vector<ComponentT*> components_;
      IdxT                    raw_nnz_{0};
    };
  }
}

#include <GridKit/Model/PhasorDynamics/ContainerImpl.hpp>
```

Copy and move are deleted because the stored pointers may refer to members of the concrete derived composite.

### `ContainerImpl.hpp`

```cpp
#pragma once

#include <algorithm>
#include <stdexcept>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/Container.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    index_type Container<scalar_type, index_type>::size()
    {
      IdxT size = 0;

      for (const auto* component : components_)
      {
        size += component->size();
      }

      return size;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::setGridKitComponentID(
        IdxT component_id)
    {
      this->gridkit_component_id_ = component_id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::addComponent(
        ComponentT* component)
    {
      if (component == nullptr)
      {
        throw std::logic_error("Container cannot contain a null component");
      }

      if (component == this)
      {
        throw std::logic_error("Container cannot contain itself");
      }

      if (this->allocated_ || this->bound_)
      {
        throw std::logic_error(
            "Container topology cannot change after allocation");
      }

      if (std::find(components_.begin(), components_.end(), component)
          != components_.end())
      {
        throw std::logic_error(
            "Container cannot contain the same component twice");
      }

      const IdxT component_id =
          static_cast<IdxT>(components_.size());

      component->setGridKitComponentID(component_id);
      component->setSystemBase(
          this->freq_system_base_,
          this->va_system_base_);

      components_.push_back(component);
    }

    template <typename scalar_type, typename index_type>
    auto Container<scalar_type, index_type>::getComponent(
        IdxT component_id) -> ComponentT*
    {
      return components_.at(static_cast<size_t>(component_id));
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::allocate()
    {
      const bool root = !this->bound_;

      this->size_ = size();

      if (root && !this->allocated_)
      {
        this->allocateVectors(this->size_);
      }

      if (this->y_.getSize() != this->size_
          || this->yp_.getSize() != this->size_
          || this->f_.getSize() != this->size_
          || this->abs_tol_.getSize() != this->size_)
      {
        throw std::runtime_error(
            "Container vector sizes do not match the container size");
      }

      this->tag_.resize(static_cast<size_t>(this->size_));
      this->variable_indices_.resize(
          static_cast<size_t>(this->size_));
      this->residual_indices_.resize(
          static_cast<size_t>(this->size_));

      IdxT offset = 0;

      for (auto* component : components_)
      {
        const IdxT child_size = component->size();

        if (component->bind(
                this->y_,
                this->yp_,
                this->f_,
                this->abs_tol_,
                offset)
            != 0)
        {
          throw std::runtime_error(
              "Container failed to bind child vectors");
        }

        if (component->allocate() != 0)
        {
          throw std::runtime_error(
              "Container failed to allocate a child component");
        }

        if (component->size() != child_size)
        {
          throw std::runtime_error(
              "A component changed size during Container allocation");
        }

        offset += child_size;
      }

      if (offset != this->size_)
      {
        throw std::runtime_error(
            "Container child sizes do not match the container size");
      }

      if (root)
      {
        for (IdxT i = 0; i < this->size_; ++i)
        {
          if (setVariableIndex(i, i) != 0
              || setResidualIndex(i, i) != 0)
          {
            throw std::runtime_error(
                "Container failed to assign root indices");
          }
        }

        if (bindSignalNodes(this->size_) != 0)
        {
          throw std::runtime_error(
              "Container failed to bind signal nodes");
        }

        const int error_count = verify();

        if (error_count != 0)
        {
          Log::error() << "Component errors: "
                       << error_count << "\n";
          throw std::runtime_error(
              "Container verification failed");
        }

        if (hasJacobian())
        {
          updateTime(RealT{0}, RealT{1});
          initialize();
          evaluateResidual();
          evaluateJacobian();
        }
      }

      this->allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::verify() const
    {
      int error_count = 0;

      for (const auto* component : components_)
      {
        error_count += component->verify();
      }

      return error_count;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::bindSignalNodes(
        IdxT system_size)
    {
      if (ComponentT::bindSignalNodes(system_size) != 0)
      {
        return 1;
      }

      for (auto* component : components_)
      {
        if (component->bindSignalNodes(system_size) != 0)
        {
          return 1;
        }
      }

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::initialize()
    {
      for (auto* component : components_)
      {
        component->initialize();
      }

      this->y_.setDataUpdated();
      this->yp_.setDataUpdated();

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::tagDifferentiable()
    {
      IdxT offset = 0;

      for (auto* component : components_)
      {
        component->tagDifferentiable();

        for (IdxT i = 0; i < component->size(); ++i)
        {
          this->tag_[static_cast<size_t>(offset + i)] =
              component->tag()[static_cast<size_t>(i)];
        }

        offset += component->size();
      }

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::setAbsoluteTolerance(
        RealT rel_tol)
    {
      for (auto* component : components_)
      {
        component->setAbsoluteTolerance(rel_tol);
      }

      this->abs_tol_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      evaluateExternalResidual();

      this->f_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::evaluateInternalResidual()
    {
      for (auto* component : components_)
      {
        component->evaluateInternalResidual();
      }

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::evaluateExternalResidual()
    {
      for (auto* component : components_)
      {
        component->evaluateExternalResidual();
      }

      return 0;
    }

    template <typename scalar_type, typename index_type>
    bool Container<scalar_type, index_type>::hasJacobian()
    {
#ifndef GRIDKIT_ENABLE_ENZYME
      return false;
#else
      for (auto* component : components_)
      {
        if (!component->hasJacobian())
        {
          return false;
        }
      }

      return true;
#endif
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::updateTime(
        RealT t,
        RealT alpha)
    {
      ComponentT::updateTime(t, alpha);

      for (auto* component : components_)
      {
        component->updateTime(t, alpha);
      }
    }

    template <typename scalar_type, typename index_type>
    void Container<scalar_type, index_type>::setSystemBase(
        RealT freq_system_base,
        RealT va_system_base)
    {
      ComponentT::setSystemBase(
          freq_system_base,
          va_system_base);

      for (auto* component : components_)
      {
        component->setSystemBase(
            freq_system_base,
            va_system_base);
      }
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::setVariableIndex(
        IdxT local_index,
        IdxT global_index)
    {
      if (ComponentT::setVariableIndex(
              local_index,
              global_index)
          != 0)
      {
        return 1;
      }

      IdxT offset = 0;

      for (auto* component : components_)
      {
        const IdxT child_size = component->size();

        if (local_index < offset + child_size)
        {
          return component->setVariableIndex(
              local_index - offset,
              global_index);
        }

        offset += child_size;
      }

      return 1;
    }

    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::setResidualIndex(
        IdxT local_index,
        IdxT global_index)
    {
      if (ComponentT::setResidualIndex(
              local_index,
              global_index)
          != 0)
      {
        return 1;
      }

      IdxT offset = 0;

      for (auto* component : components_)
      {
        const IdxT child_size = component->size();

        if (local_index < offset + child_size)
        {
          return component->setResidualIndex(
              local_index - offset,
              global_index);
        }

        offset += child_size;
      }

      return 1;
    }
```

The Jacobian implementation preserves a raw global-index COO for parents and builds a deduplicated CSR only at the direct root:

```cpp
    template <typename scalar_type, typename index_type>
    int Container<scalar_type, index_type>::evaluateJacobian()
    {
      IdxT new_raw_nnz = 0;

      for (auto* component : components_)
      {
        component->evaluateJacobian();

        auto* jacobian = component->getCooJacobian();

        if (jacobian == nullptr)
        {
          if (component->nnz() != IdxT{})
          {
            throw std::runtime_error(
                "Container child returned a null COO Jacobian");
          }

          continue;
        }

        new_raw_nnz += jacobian->getNnz();
      }

      if (this->coo_jac_ == nullptr)
      {
        raw_nnz_ = new_raw_nnz;
        this->nnz_ = raw_nnz_;

        this->J_rows_buffer_ =
            new IdxT[static_cast<size_t>(raw_nnz_)];
        this->J_cols_buffer_ =
            new IdxT[static_cast<size_t>(raw_nnz_)];
        this->J_vals_buffer_ =
            new RealT[static_cast<size_t>(raw_nnz_)];

        IdxT counter  = 0;
        IdxT num_rows = 0;
        IdxT num_cols = 0;

        for (auto* component : components_)
        {
          auto* jacobian = component->getCooJacobian();

          if (jacobian == nullptr)
          {
            continue;
          }

          const auto* rows = jacobian->getRowData();
          const auto* cols = jacobian->getColData();
          const auto* vals = jacobian->getValues();

          for (IdxT i = 0; i < jacobian->getNnz(); ++i)
          {
            if (rows[i] == INVALID_INDEX<IdxT>
                || cols[i] == INVALID_INDEX<IdxT>)
            {
              throw std::runtime_error(
                  "Container child returned an invalid Jacobian index");
            }

            this->J_rows_buffer_[counter] = rows[i];
            this->J_cols_buffer_[counter] = cols[i];
            this->J_vals_buffer_[counter] = vals[i];

            if (rows[i] + 1 > num_rows)
            {
              num_rows = rows[i] + 1;
            }

            if (cols[i] + 1 > num_cols)
            {
              num_cols = cols[i] + 1;
            }

            ++counter;
          }
        }

        this->coo_jac_ =
            new CooMatrixT(num_rows, num_cols, raw_nnz_);

        if (this->coo_jac_->setDataPointers(
                this->J_rows_buffer_,
                this->J_cols_buffer_,
                this->J_vals_buffer_,
                memory::HOST)
            != 0)
        {
          throw std::runtime_error(
              "Container failed to construct its COO Jacobian");
        }
      }
      else
      {
        if (new_raw_nnz != raw_nnz_)
        {
          throw std::runtime_error(
              "Container Jacobian nonzero count changed after allocation");
        }

        IdxT counter = 0;

        for (auto* component : components_)
        {
          auto* jacobian = component->getCooJacobian();

          if (jacobian == nullptr)
          {
            continue;
          }

          const auto* rows = jacobian->getRowData();
          const auto* cols = jacobian->getColData();
          const auto* vals = jacobian->getValues();

          for (IdxT i = 0; i < jacobian->getNnz(); ++i)
          {
            if (this->J_rows_buffer_[counter] != rows[i]
                || this->J_cols_buffer_[counter] != cols[i])
            {
              throw std::runtime_error(
                  "Container Jacobian pattern changed after allocation");
            }

            this->J_vals_buffer_[counter] = vals[i];
            ++counter;
          }
        }

        this->nnz_ = raw_nnz_;
      }

      // A nested Container exports its unsorted, undeduplicated global COO.
      if (this->bound_)
      {
        return 0;
      }

      // A direct root must be closed.
      for (IdxT i = 0; i < raw_nnz_; ++i)
      {
        if (this->J_rows_buffer_[i] >= this->size_
            || this->J_cols_buffer_[i] >= this->size_)
        {
          throw std::runtime_error(
              "A root Container Jacobian references an external variable");
        }
      }

      if (this->csr_jac_ == nullptr)
      {
        auto* rows =
            new IdxT[static_cast<size_t>(raw_nnz_)];
        auto* cols =
            new IdxT[static_cast<size_t>(raw_nnz_)];
        auto* vals =
            new RealT[static_cast<size_t>(raw_nnz_)];

        for (IdxT i = 0; i < raw_nnz_; ++i)
        {
          rows[i] = this->J_rows_buffer_[i];
          cols[i] = this->J_cols_buffer_[i];
          vals[i] = this->J_vals_buffer_[i];
        }

        CooMatrixT jacobian(
            this->size_,
            this->size_,
            raw_nnz_,
            &rows,
            &cols,
            &vals);

        IdxT* row_ptrs = jacobian.getCsrRowData();

        this->nnz_ = jacobian.getNnz();

        auto* csr_cols =
            new IdxT[static_cast<size_t>(this->nnz_)];
        auto* csr_vals =
            new RealT[static_cast<size_t>(this->nnz_)];

        for (IdxT i = 0; i < this->nnz_; ++i)
        {
          csr_cols[i] = jacobian.getColData()[i];
          csr_vals[i] = jacobian.getValues()[i];
        }

        this->csr_jac_ = new typename ComponentT::CsrMatrixT(
            this->size_,
            this->size_,
            this->nnz_,
            &row_ptrs,
            &csr_cols,
            &csr_vals);

        const auto* map_to_sorted =
            jacobian.getMapToSorted();
        const auto* map_to_deduplicated =
            jacobian.getMapToDeduplicated();

        this->map_to_csr_ =
            new IdxT[static_cast<size_t>(raw_nnz_)];

        for (IdxT i = 0; i < raw_nnz_; ++i)
        {
          this->map_to_csr_[map_to_sorted[i]] =
              map_to_deduplicated[i];
        }
      }
      else
      {
        auto* values = this->csr_jac_->getValues();

        for (IdxT i = 0; i < this->csr_jac_->getNnz(); ++i)
        {
          values[i] = RealT{0};
        }

        for (IdxT i = 0; i < raw_nnz_; ++i)
        {
          values[this->map_to_csr_[i]] +=
              this->J_vals_buffer_[i];
        }

        this->csr_jac_->setUpdated(memory::HOST);
      }

      return 0;
    }
  }
}
```

Child COO indices are already global and are copied unchanged. The root CSR is built from a temporary COO because `getCsrRowData()` sorts and deduplicates its input in place; calling it on the retained aggregate COO would break nesting.

### Core CMake registration

Add only the two headers:

```cmake
gridkit_add_library(
  phasor_dynamics_core
  INTERFACE_TARGET
  HEADERS BusBase.hpp
          Component.hpp
          ComponentData.hpp
          ComponentLibrary.hpp
          ComponentSignals.hpp
          Container.hpp
          ContainerImpl.hpp
  # existing link libraries remain unchanged
)
```

Do not add `Container` to `ComponentLibrary.hpp`; callers include `Container.hpp` directly.

## 4. Container test model and acceptance tests

### Test component

The test component supplies one owned variable \(x\), reads external \(u\), and implements:

\[
F_x = \dot{x} + a x - b u,\qquad F_u \mathrel{+}= c x.
\]

It splits the state and derivative diagonal stamps to exercise root COO deduplication.

```cpp
enum class LinearInternalVariables : size_t
{
  X,
  MAXIMUM
};

enum class LinearExternalVariables : size_t
{
  U,
  MAXIMUM
};

template <typename ScalarT, typename IdxT>
class LinearComponent
  : public PhasorDynamics::Component<ScalarT, IdxT>
{
public:
  using ComponentT = PhasorDynamics::Component<ScalarT, IdxT>;
  using RealT      = typename ComponentT::RealT;
  using SignalsT   = PhasorDynamics::ComponentSignals<
      ScalarT,
      IdxT,
      LinearInternalVariables,
      LinearExternalVariables>;

  LinearComponent(
      RealT   a,
      RealT   b,
      RealT   c,
      ScalarT y0,
      ScalarT yp0)
    : a_(a),
      b_(b),
      c_(c),
      y0_(y0),
      yp0_(yp0)
  {
    this->size_ = 1;
    this->nnz_  = 4;
  }

  int setGridKitComponentID(IdxT component_id) override
  {
    this->gridkit_component_id_ = component_id;
    return 0;
  }

  int allocate() override
  {
    if (!this->allocated_)
    {
      this->allocateVectors(this->size_);
    }

    this->tag_.resize(1);
    this->variable_indices_.resize(1);
    this->residual_indices_.resize(1);

    this->setVariableIndex(0, 0);
    this->setResidualIndex(0, 0);

    if (signals_.template isAssigned<LinearInternalVariables::X>())
    {
      signals_.template getSignalNode<LinearInternalVariables::X>()
          ->set(
              this->y_int_,
              this->yp_int_,
              this->f_int_,
              &this->getVariableIndex(0),
              &this->getResidualIndex(0));
    }

    if (signals_.template isAttached<LinearExternalVariables::U>())
    {
      if (ComponentT::bind(
              signals_.template getSignalNode<
                  LinearExternalVariables::U>(),
              0)
          != 0)
      {
        return 1;
      }
    }

    this->allocated_ = true;
    return 0;
  }

  int verify() const override
  {
    int errors = 0;

    if (!signals_.template isAssigned<LinearInternalVariables::X>()
        || !signals_.template getSignalNode<
                LinearInternalVariables::X>()->linked())
    {
      ++errors;
    }

    if (!signals_.template isAttached<LinearExternalVariables::U>()
        || !signals_.template isLinked<LinearExternalVariables::U>())
    {
      ++errors;
    }

    return errors;
  }

  int initialize() override
  {
    this->y_int_[0]  = y0_;
    this->yp_int_[0] = yp0_;
    return 0;
  }

  int tagDifferentiable() override
  {
    this->tag_[0] = true;
    return 0;
  }

  int setAbsoluteTolerance(RealT rel_tol) override
  {
    this->abs_tol_.getData()[0] = rel_tol;
    this->abs_tol_.setDataUpdated();
    return 0;
  }

  int evaluateResidual() override
  {
    evaluateInternalResidual();
    evaluateExternalResidual();
    this->f_.setDataUpdated();
    return 0;
  }

  int evaluateInternalResidual() override
  {
    this->f_int_[0] =
        this->yp_int_[0]
        + a_ * this->y_int_[0]
        - b_ * *this->y_ext_[0];

    return 0;
  }

  int evaluateExternalResidual() override
  {
    *this->f_ext_[0] += c_ * this->y_int_[0];
    return 0;
  }

  int evaluateJacobian() override
  {
    if (this->J_rows_buffer_ == nullptr)
    {
      this->J_rows_buffer_ = new IdxT[4];
      this->J_cols_buffer_ = new IdxT[4];
      this->J_vals_buffer_ = new RealT[4];
    }

    const IdxT x_variable = this->getVariableIndex(0);
    const IdxT x_residual = this->getResidualIndex(0);
    const IdxT u_variable = *this->variable_indices_ext_[0];
    const IdxT u_residual = *this->residual_indices_ext_[0];

    this->J_rows_buffer_[0] = x_residual;
    this->J_cols_buffer_[0] = x_variable;
    this->J_vals_buffer_[0] = a_;

    this->J_rows_buffer_[1] = x_residual;
    this->J_cols_buffer_[1] = x_variable;
    this->J_vals_buffer_[1] = this->alpha_;

    this->J_rows_buffer_[2] = x_residual;
    this->J_cols_buffer_[2] = u_variable;
    this->J_vals_buffer_[2] = -b_;

    this->J_rows_buffer_[3] = u_residual;
    this->J_cols_buffer_[3] = x_variable;
    this->J_vals_buffer_[3] = c_;

    this->constructCoo();
    return 0;
  }

  auto getSignals() -> SignalsT&
  {
    return signals_;
  }

private:
  SignalsT signals_;

  RealT   a_;
  RealT   b_;
  RealT   c_;
  ScalarT y0_;
  ScalarT yp0_;
};
```

The external derivative rule is explicit: if an equation reads `yp_ext_`, that component stamps `alpha * dF/dyp_ext` itself. `SignalNode` carries aliases but owns no equation or Jacobian behavior.

### Dedicated tests

Create `ContainerTests.hpp` and `runContainerTests.cpp`. The main composition test uses:

```cpp
TestOutcome composition()
{
  TestStatus success = true;

  ClosedModel<ScalarT, IdxT> model;
  Model::Evaluator<ScalarT, IdxT>* evaluator = &model;

  success *= evaluator->allocate() == 0;
  success *= evaluator->size() == 3;
  success *= evaluator->initialize() == 0;

  evaluator->updateTime(RealT{0}, RealT{10});

  success *= evaluator->tagDifferentiable() == 0;
  success *= evaluator->setAbsoluteTolerance(RealT{1.0e-6}) == 0;
  success *= evaluator->evaluateResidual() == 0;

  const ScalarT expected_residual[] = {
      ScalarT{-27.9},
      ScalarT{41.2},
      ScalarT{13.3}};

  const auto* residual = evaluator->getResidual().getData();

  for (IdxT i = 0; i < 3; ++i)
  {
    success *= isEqual(residual[i], expected_residual[i]);
  }

  auto& cascade = model.cascade();

  success *= cascade.getSignals()
                 .template getSignalNode<CascadeExternalVariables::U>()
             == model.uNode();

  success *= cascade.getSignals()
                 .template getSignalNode<CascadeInternalVariables::X>()
             == model.zNode();

  success *= cascade.y().getData()
             == model.y().getData() + 1;

  success *= cascade.getComponent(0)->y().getData()
             == model.y().getData() + 1;

  success *= cascade.getComponent(1)->y().getData()
             == model.y().getData() + 2;

  success *= evaluator->evaluateJacobian() == 0;

  const auto* raw = model.getCooJacobian();
  const auto* csr = evaluator->getCsrJacobian();

  success *= raw != nullptr;
  success *= raw->getNnz() == 12;
  success *= csr != nullptr;
  success *= csr->getNumRows() == 3;
  success *= csr->getNumColumns() == 3;
  success *= csr->getNnz() == 9;

  const IdxT expected_rows[] = {0, 3, 6, 9};
  const IdxT expected_cols[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  const RealT expected_values[] = {
      11.0, 5.0, -9.0,
      -4.0, 13.0, 8.0,
      2.0, -7.0, 16.0};

  for (IdxT i = 0; i < 4; ++i)
  {
    success *= csr->getRowData()[i] == expected_rows[i];
  }

  for (IdxT i = 0; i < 9; ++i)
  {
    success *= csr->getColData()[i] == expected_cols[i];
    success *= isEqual(csr->getValues()[i], expected_values[i]);
  }

  return success.report(__func__);
}
```

This expected residual distinguishes correct two-phase evaluation from child-by-child assignment/addition.

Test nesting at a nonzero `SystemModel` offset:

```cpp
TestOutcome nestedInSystemModel()
{
  TestStatus success = true;

  PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
  ClosedModel<ScalarT, IdxT>         model;
  PhasorDynamics::SystemModel<ScalarT, IdxT> system;

  system.addBus(&bus);
  system.addComponent(&model);

  success *= system.allocate() == 0;
  success *= system.initialize() == 0;

  system.updateTime(RealT{0}, RealT{10});

  success *= system.evaluateResidual() == 0;

  const ScalarT expected_residual[] = {
      ScalarT{-27.9},
      ScalarT{41.2},
      ScalarT{13.3}};

  for (IdxT i = 0; i < 3; ++i)
  {
    success *= isEqual(
        model.getResidual().getData()[i],
        expected_residual[i]);
  }

  const IdxT offset = bus.size();

  for (IdxT i = 0; i < model.size(); ++i)
  {
    success *= model.getVariableIndex(i) == offset + i;
    success *= model.getResidualIndex(i) == offset + i;
  }

  success *= system.evaluateJacobian() == 0;

  const auto* raw = model.getCooJacobian();
  success *= raw != nullptr;
  success *= raw->getNnz() == 12;

  for (IdxT i = 0; i < raw->getNnz(); ++i)
  {
    success *= raw->getRowData()[i] >= offset;
    success *= raw->getRowData()[i] < offset + model.size();
    success *= raw->getColData()[i] >= offset;
    success *= raw->getColData()[i] < offset + model.size();
  }

  return success.report(__func__);
}
```

Test zero-state composition and topology locking:

```cpp
TestOutcome zeroStateAndTopology()
{
  TestStatus success = true;

  PhasorDynamics::Container<ScalarT, IdxT> child;
  PhasorDynamics::Container<ScalarT, IdxT> root;
  PhasorDynamics::Container<ScalarT, IdxT> late;

  root.addComponent(&child);

  success *= root.allocate() == 0;
  success *= root.size() == 0;
  success *= child.size() == 0;
  success *= root.initialize() == 0;
  success *= root.evaluateResidual() == 0;
  success *= root.evaluateJacobian() == 0;

  bool late_add_threw = false;

  try
  {
    root.addComponent(&late);
  }
  catch (const std::logic_error&)
  {
    late_add_threw = true;
  }

  success *= late_add_threw;

  return success.report(__func__);
}
```

Test that a legacy node remains usable as an ordinary signal but cannot populate the complete external arrays:

```cpp
TestOutcome incompleteExternalBinding()
{
  TestStatus success = true;

  ScalarT external_value{1.0};
  IdxT    external_index{0};

  PhasorDynamics::SignalNode<ScalarT, IdxT> input;
  PhasorDynamics::SignalNode<ScalarT, IdxT> output;

  input.set(&external_value, &external_index);

  LinearComponent<ScalarT, IdxT> component(
      1.0, 1.0, 1.0, 0.0, 0.0);

  component.getSignals()
      .template attachSignalNode<LinearExternalVariables::U>(&input);
  component.getSignals()
      .template assignSignalNode<LinearInternalVariables::X>(&output);

  PhasorDynamics::SystemModel<ScalarT, IdxT> system;
  system.addComponent(&component);

  bool allocation_threw = false;

  try
  {
    system.allocate();
  }
  catch (const std::runtime_error&)
  {
    allocation_threw = true;
  }

  success *= allocation_threw;
  success *= input.linked();
  success *= isEqual(input.read(), external_value);

  return success.report(__func__);
}
```

Runner:

```cpp
#include "ContainerTests.hpp"

int main()
{
  using namespace GridKit::Testing;

  TestingResults result;
  ContainerTests<double, size_t> test;

  result += test.composition();
  result += test.nestedInSystemModel();
  result += test.zeroStateAndTopology();
  result += test.incompleteExternalBinding();

  return result.summary();
}
```

CMake:

```cmake
add_executable(test_phasor_container runContainerTests.cpp)
target_link_libraries(
  test_phasor_container
  GridKit::definitions
  GridKit::phasor_dynamics_systemmodel
  GridKit::testing)

add_test(
  NAME PhasorDynamicsContainerTest
  COMMAND test_phasor_container)
```

Add `test_phasor_container` to the existing unit-test install target list.

Verification:

```bash
cmake --build build --target test_phasor_system -j 10
cmake --build build --target test_phasor_system_single_component -j 10
cmake --build build --target test_phasor_container -j 10
ctest --test-dir build -R 'PhasorDynamics(System|Container)' --output-on-failure
cmake --build build --target DynamicSimulation -j 10
```

Acceptance requires:

- Existing two-pointer `SignalNode` users compile and behave unchanged.
- Complete nodes alias `y`, `yp`, `f`, and both global indices.
- Cyclic sibling connections resolve after every producer has allocated.
- Every internal residual assignment precedes every external `+=`.
- Public ports use the exact caller-supplied nodes.
- Root, nested, child, and zero-state vectors alias correctly.
- Nested COO indices remain global and are not offset twice.
- Root COO remains raw; root CSR is sorted and deduplicated.
- A direct root rejects Jacobian references outside its owned variables.
- Both current sparse/dense solver-selection configurations compile.
- Existing PhasorDynamics system tests and `DynamicSimulation` still build.

## Assumptions and explicit exclusions

- This targets the current `ComponentSignals` API. It does not depend on or reproduce the open IOPorts work.
- The `*_int_`/`*_ext_` split adopts the useful ownership semantics of the open external-variable work, but does not copy its model-specific implementation.
- Existing component output calls remain on `SignalNode::set(y, index)` in this change.
- Only new components that require derivative/residual aliases use the five-pointer overload and overloaded `Component::bind`.
- `SystemModel` is not changed to inherit `Container`; that becomes a small later migration after this design is established.
- A bus-coupled component can be contained when the `Container` is nested under `SystemModel`; it cannot make a standalone root closed because `Container` does not contain buses.
- Child traversal keeps current `SystemModel` status behavior: `verify()` accumulates errors, allocation checks structural failures, and ordinary initialize/residual/Jacobian traversal returns zero.
- Signal rewiring after allocation is unsupported by contract; no new locking state is added to `ComponentSignals`.
- No changes are made to `Vector`, buses, monitors, parsers, model factories, existing component implementations, or `ComponentLibrary.hpp`.
