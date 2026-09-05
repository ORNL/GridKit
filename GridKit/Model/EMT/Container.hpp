#pragma once

#include <concepts>
#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <GridKit/Model/EMT/Component.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct ContainerData;

    /**
     * @brief A component that owns and composes child components.
     *
     * A Container owns no variables or equations of its own. Its state is the
     * contiguous concatenation of its children's state, and its residual and
     * Jacobian traversals preserve the internal-then-external assembly phases
     * across the complete subtree.
     */
    template <typename scalar_type, typename index_type>
    class Container : public Component<scalar_type, index_type>
    {
    protected:
      using Component<scalar_type, index_type>::abs_tol_;
      using Component<scalar_type, index_type>::allocated_;
      using Component<scalar_type, index_type>::alpha_;
      using Component<scalar_type, index_type>::coo_jac_;
      using Component<scalar_type, index_type>::csr_jac_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::map_to_csr_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::time_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using CsrMatrixT = typename Component<ScalarT, IdxT>::CsrMatrixT;
      using CooMatrixT = typename Component<ScalarT, IdxT>::CooMatrixT;
      using VectorT    = typename Component<ScalarT, IdxT>::VectorT;
      using ComponentT = Component<ScalarT, IdxT>;
      using SignalT    = Signal<ScalarT, IdxT>;
      using Port3T     = Port3<ScalarT, IdxT>;
      using ModelDataT = ContainerData<RealT, IdxT>;

      Container();
      explicit Container(const ModelDataT& data);
      ~Container() override = default;

      /** Add one owned child to this lexical scope. */
      ComponentT& add(std::string id, std::unique_ptr<ComponentT> child);

      template <typename ChildT, typename... Args>
        requires std::derived_from<ChildT, ComponentT>
      ChildT& add(std::string id, Args&&... args)
      {
        auto  child = std::make_unique<ChildT>(std::forward<Args>(args)...);
        auto* ptr   = child.get();
        add(std::move(id), std::move(child));
        return *ptr;
      }

      /** Add one owned scalar signal to this lexical scope. */
      SignalT& addSignal(std::string id);

      ComponentT&       component(std::string_view path);
      const ComponentT& component(std::string_view path) const;

      template <typename ChildT>
      ChildT& component(std::string_view path)
      {
        auto* child = dynamic_cast<ChildT*>(&component(path));
        if (child == nullptr)
        {
          throw std::invalid_argument("Component \"" + std::string(path)
                                      + "\" has the wrong class");
        }
        return *child;
      }

      template <typename ChildT>
      const ChildT& component(std::string_view path) const
      {
        auto* child = dynamic_cast<const ChildT*>(&component(path));
        if (child == nullptr)
        {
          throw std::invalid_argument("Component \"" + std::string(path)
                                      + "\" has the wrong class");
        }
        return *child;
      }

      ComponentT&       component(IdxT local_index);
      const ComponentT& component(IdxT local_index) const;

      SignalT&       signal(std::string_view path);
      const SignalT& signal(std::string_view path) const;

      /** Bind a public input to its parent-scope endpoint. */
      void input(std::string name, SignalT& signal);
      void input(std::string name, Port3T& port);

      SignalT&       inputSignal(std::string_view name);
      const SignalT& inputSignal(std::string_view name) const;
      Port3T&        inputPort(std::string_view name);
      const Port3T&  inputPort(std::string_view name) const;

      /** Define a public scalar or three-phase output. */
      void output(std::string name, SignalT& signal);
      void output(std::string name, Port3T& port);

      SignalT&       outputSignal(std::string_view name);
      const SignalT& outputSignal(std::string_view name) const;
      Port3T&        outputPort(std::string_view name);
      const Port3T&  outputPort(std::string_view name) const;

      IdxT size() override;
      int  bind(VectorT& y,
                VectorT& yp,
                VectorT& f,
                VectorT& abs_tol,
                IdxT     offset) override;
      int  assignGlobalIndices(IdxT first) override;

      int  setGridKitComponentID(IdxT component_id) override;
      int  allocate() override;
      int  verify() const override;
      int  initialize() override;
      int  tagDifferentiable() override;
      int  setAbsoluteTolerance(RealT rel_tol) override;
      int  evaluateInternalResidual() override;
      int  evaluateExternalResidual() override;
      int  evaluateResidual() override;
      int  evaluateJacobian() override;
      bool hasJacobian() override;

      void updateTime(RealT t, RealT a) override;
      void resetJacobianStructure() override;

    protected:
      bool boundToParent() const noexcept
      {
        return bound_;
      }

      template <typename Visitor>
      void forEachComponent(Visitor&& visitor)
      {
        for (auto& child : children_)
        {
          visitor(*child);
          if (auto* container = dynamic_cast<Container*>(child.get()))
          {
            container->forEachComponent(visitor);
          }
        }
      }

      template <typename Visitor>
      void forEachComponent(Visitor&& visitor) const
      {
        for (const auto& child : children_)
        {
          visitor(*child);
          if (const auto* container = dynamic_cast<const Container*>(child.get()))
          {
            container->forEachComponent(visitor);
          }
        }
      }

    private:
      using Endpoint = std::variant<SignalT*, Port3T*>;

      static void validateName(std::string_view name, std::string_view kind);
      void        declare(const ModelDataT& data, std::string path);
      void        assemble(const ModelDataT& data);
      void        wire(const ModelDataT& data);
      void        validateBoundary() const;
      void        declareInput(std::string name);
      void        bindInput(std::string_view name, Endpoint endpoint);
      void        refreshLayout();

      std::string qualify(std::string_view local_name) const;

      SignalT& source(std::string_view reference);
      Port3T&  port(std::string_view reference);
      Endpoint endpoint(std::string_view reference);
      Endpoint resolveOutput(std::string_view reference);
      Endpoint inputEndpoint(std::string_view name) const;
      Endpoint outputEndpoint(std::string_view name);

      Container&       childContainer(std::string_view id);
      const Container& childContainer(std::string_view id) const;

      // Children are destroyed before the Signals to which they may refer.
      std::vector<std::unique_ptr<SignalT>>    signals_;
      std::vector<std::unique_ptr<ComponentT>> children_;

      std::map<std::string, SignalT*, std::less<>>    signals_by_id_;
      std::map<std::string, ComponentT*, std::less<>> children_by_id_;

      std::set<std::string, std::less<>>           input_names_;
      std::map<std::string, Endpoint, std::less<>> inputs_;
      std::map<std::string, Endpoint, std::less<>> outputs_;

      std::vector<IdxT> offsets_;
      std::string       path_;
      size_t            jacobian_capacity_{0};
      bool              bound_{false};
    };
  } // namespace EMT
} // namespace GridKit
