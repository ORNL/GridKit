/**
 * @file ComponentLibrary.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT models that a case can name.
 */

#pragma once

#include <type_traits>

#include <GridKit/Model/EMT/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Component/Controller/Esdc1a/Esdc1a.hpp>
#include <GridKit/Model/EMT/Component/Controller/Esst1a/Esst1a.hpp>
#include <GridKit/Model/EMT/Component/Controller/Esst4b/Esst4b.hpp>
#include <GridKit/Model/EMT/Component/Controller/Exst1Pti/Exst1Pti.hpp>
#include <GridKit/Model/EMT/Component/Controller/GastPti/GastPti.hpp>
#include <GridKit/Model/EMT/Component/Controller/Ggov1/Ggov1.hpp>
#include <GridKit/Model/EMT/Component/Controller/Hygov/Hygov.hpp>
#include <GridKit/Model/EMT/Component/Controller/Ieeeg1/Ieeeg1.hpp>
#include <GridKit/Model/EMT/Component/Controller/Ieeest/Ieeest.hpp>
#include <GridKit/Model/EMT/Component/Controller/Ieeet1/Ieeet1.hpp>
#include <GridKit/Model/EMT/Component/Controller/Reeca/Reeca.hpp>
#include <GridKit/Model/EMT/Component/Controller/Reecb/Reecb.hpp>
#include <GridKit/Model/EMT/Component/Controller/Repca/Repca.hpp>
#include <GridKit/Model/EMT/Component/Controller/SexsPti/SexsPti.hpp>
#include <GridKit/Model/EMT/Component/Controller/Tgov1/Tgov1.hpp>
#include <GridKit/Model/EMT/Component/Converter/Regca/Regca.hpp>
#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumped.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadRLC/LoadRLC.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/EMT/Component/Source/Machine/Machine.hpp>
#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSource.hpp>
#include <GridKit/Model/EMT/Component/Switch/Switch.hpp>
#include <GridKit/Model/EMT/Component/Transformer/Transformer.hpp>
#include <GridKit/Model/EMT/Component/TransformerThreeWinding/TransformerThreeWinding.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Park/Park.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Pll/Pll.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Pmu/Pmu.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Signal/Signal.hpp>

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
        visit(std::type_identity<Esdc1a<ScalarT, IdxT>>{});
        visit(std::type_identity<Esst1a<ScalarT, IdxT>>{});
        visit(std::type_identity<Esst4b<ScalarT, IdxT>>{});
        visit(std::type_identity<Exst1Pti<ScalarT, IdxT>>{});
        visit(std::type_identity<GastPti<ScalarT, IdxT>>{});
        visit(std::type_identity<Ggov1<ScalarT, IdxT>>{});
        visit(std::type_identity<Hygov<ScalarT, IdxT>>{});
        visit(std::type_identity<Ieeeg1<ScalarT, IdxT>>{});
        visit(std::type_identity<Ieeest<ScalarT, IdxT>>{});
        visit(std::type_identity<Ieeet1<ScalarT, IdxT>>{});
        visit(std::type_identity<LineLumped<ScalarT, IdxT>>{});
        visit(std::type_identity<LoadRLC<ScalarT, IdxT>>{});
        visit(std::type_identity<LoadZ<ScalarT, IdxT>>{});
        visit(std::type_identity<Machine<ScalarT, IdxT>>{});
        visit(std::type_identity<Park<ScalarT, IdxT>>{});
        visit(std::type_identity<Pll<ScalarT, IdxT>>{});
        visit(std::type_identity<Pmu<ScalarT, IdxT>>{});
        visit(std::type_identity<Reeca<ScalarT, IdxT>>{});
        visit(std::type_identity<Reecb<ScalarT, IdxT>>{});
        visit(std::type_identity<Regca<ScalarT, IdxT>>{});
        visit(std::type_identity<Repca<ScalarT, IdxT>>{});
        visit(std::type_identity<SexsPti<ScalarT, IdxT>>{});
        visit(std::type_identity<Signal<ScalarT, IdxT>>{});
        visit(std::type_identity<Switch<ScalarT, IdxT>>{});
        visit(std::type_identity<Tgov1<ScalarT, IdxT>>{});
        visit(std::type_identity<Transformer<ScalarT, IdxT>>{});
        visit(std::type_identity<TransformerThreeWinding<ScalarT, IdxT>>{});
        visit(std::type_identity<VectorFit<ScalarT, IdxT>>{});
        visit(std::type_identity<VoltageSource<ScalarT, IdxT>>{});
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
