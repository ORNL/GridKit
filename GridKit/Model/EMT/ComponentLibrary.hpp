/**
 * @file ComponentLibrary.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT models that a case can name.
 */

#pragma once

#include <type_traits>

#include <GridKit/Model/EMT/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZ.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename ScalarT, typename IdxT>
    struct ComponentLibrary
    {
      template <typename VisitorT>
      static void forEach(VisitorT&& visit)
      {
        visit(std::type_identity<Bus<ScalarT, IdxT>>{});
        visit(std::type_identity<LoadZ<ScalarT, IdxT>>{});
      }

      /// Calls `visitor` with `component` as its library model
      template <typename VisitorT>
      static void visit(Component<ScalarT, IdxT>& component, VisitorT&& visitor)
      {
        forEach(
            [&]<typename ModelT>(std::type_identity<ModelT>)
            {
              if (auto* model = dynamic_cast<ModelT*>(&component))
              {
                visitor(*model);
              }
            });
      }
    };
  } // namespace EMT
} // namespace GridKit
