#include <cassert>
#include <iostream>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusFactory.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>

// Include all components
#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for the system model
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel()
      : monitor_(std::make_unique<MonitorT>())
    {
    }

    /**
     * @brief Construct a new System Model object
     *
     * @param[in] data - Data structure with complete system data
     *
     * @pre SystemModelData contains consistent connectivity information
     * and physically meaningful model parameters.
     *
     * @post All component models in SystemModelData are created, and
     * correctly connected into the system model.
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel(const SystemModelData<RealT, IdxT>& data)
      : monitor_(std::make_unique<MonitorT>(time_))
    {
      using namespace Governor;
      using namespace Exciter;
      using namespace Stabilizer;
      using namespace Controller;
      using namespace Converter;

      owns_components_ = true;

      // Store parsed system bases before constructing data-driven components.
      this->setSystemBase(data.freq_base, data.va_base);

      // Add electrical buses
      for (const auto& busdata : data.bus)
      {
        BusBase<ScalarT, IdxT>* bus = BusFactory<ScalarT, IdxT>::create(busdata);
        addBus(bus);
      }

      /// @todo Signal data needs to be populated by the parser.
      /// See TODOs in SystemModelDataJSONParser
      for (const auto& signaldata : data.signal)
      {
        SignalNode<ScalarT, IdxT>* signal = new SignalNode<ScalarT, IdxT>(signaldata);
        addSignal(signal);
      }

      // Add bus-to-signal adapters
      for (const auto& adapterdata : data.adapter)
      {
        IdxT bus_index = 0;
        if (adapterdata.buses.contains(BusToSignalAdapterBuses::bus))
        {
          bus_index = adapterdata.buses.at(BusToSignalAdapterBuses::bus);
        }

        auto* adapter = new BusToSignalAdapter<ScalarT, IdxT>(getBus(bus_index));

        if (adapterdata.signal_outputs.contains(BusToSignalAdapterSignalOutputs::vr))
        {
          IdxT           vr    = adapterdata.signal_outputs.at(BusToSignalAdapterSignalOutputs::vr);
          constexpr auto VREAL = BusToSignalAdapterInternalVariables::VREAL;
          adapter->getSignals().template assignSignalNode<VREAL>(getSignal(vr));
        }

        if (adapterdata.signal_outputs.contains(BusToSignalAdapterSignalOutputs::vi))
        {
          IdxT           vi    = adapterdata.signal_outputs.at(BusToSignalAdapterSignalOutputs::vi);
          constexpr auto VIMAG = BusToSignalAdapterInternalVariables::VIMAG;
          adapter->getSignals().template assignSignalNode<VIMAG>(getSignal(vi));
        }

        if (adapterdata.signal_inputs.contains(BusToSignalAdapterSignalInputs::ir))
        {
          IdxT           ir    = adapterdata.signal_inputs.at(BusToSignalAdapterSignalInputs::ir);
          constexpr auto IREAL = BusToSignalAdapterExternalVariables::IREAL;
          adapter->getSignals().template attachSignalNode<IREAL>(getSignal(ir));
        }

        if (adapterdata.signal_inputs.contains(BusToSignalAdapterSignalInputs::ii))
        {
          IdxT           ii    = adapterdata.signal_inputs.at(BusToSignalAdapterSignalInputs::ii);
          constexpr auto IIMAG = BusToSignalAdapterExternalVariables::IIMAG;
          adapter->getSignals().template attachSignalNode<IIMAG>(getSignal(ii));
        }

        addComponent(adapter);
      }

      // Add REGCA converters
      for (const auto& regcadata : data.regca)
      {
        IdxT bus_index = 0;
        if (regcadata.buses.contains(RegcaBuses::bus))
        {
          bus_index = regcadata.buses.at(RegcaBuses::bus);
        }

        auto* regca = new Regca<ScalarT, IdxT>(getBus(bus_index), regcadata);

        if (regcadata.signal_inputs.contains(RegcaSignalInputs::ipcmd))
        {
          IdxT           ipcmd = regcadata.signal_inputs.at(RegcaSignalInputs::ipcmd);
          constexpr auto IPCMD = RegcaExternalVariables::IPCMD;
          regca->getSignals().template attachSignalNode<IPCMD>(getSignal(ipcmd));
        }

        if (regcadata.signal_inputs.contains(RegcaSignalInputs::iqcmd))
        {
          IdxT           iqcmd = regcadata.signal_inputs.at(RegcaSignalInputs::iqcmd);
          constexpr auto IQCMD = RegcaExternalVariables::IQCMD;
          regca->getSignals().template attachSignalNode<IQCMD>(getSignal(iqcmd));
        }

        if (regcadata.signal_outputs.contains(RegcaSignalOutputs::ibranchr))
        {
          IdxT           ibranchr = regcadata.signal_outputs.at(RegcaSignalOutputs::ibranchr);
          constexpr auto IR       = RegcaInternalVariables::IR;
          regca->getSignals().template assignSignalNode<IR>(getSignal(ibranchr));
        }

        if (regcadata.signal_outputs.contains(RegcaSignalOutputs::ibranchi))
        {
          IdxT           ibranchi = regcadata.signal_outputs.at(RegcaSignalOutputs::ibranchi);
          constexpr auto II       = RegcaInternalVariables::II;
          regca->getSignals().template assignSignalNode<II>(getSignal(ibranchi));
        }

        if (regcadata.signal_outputs.contains(RegcaSignalOutputs::pbranch))
        {
          IdxT           pbranch = regcadata.signal_outputs.at(RegcaSignalOutputs::pbranch);
          constexpr auto PBR     = RegcaInternalVariables::PBR;
          regca->getSignals().template assignSignalNode<PBR>(getSignal(pbranch));
        }

        if (regcadata.signal_outputs.contains(RegcaSignalOutputs::qbranch))
        {
          IdxT           qbranch = regcadata.signal_outputs.at(RegcaSignalOutputs::qbranch);
          constexpr auto QBR     = RegcaInternalVariables::QBR;
          regca->getSignals().template assignSignalNode<QBR>(getSignal(qbranch));
        }

        addComponent(regca);
      }

      // Add branches
      for (const auto& branchdata : data.branch)
      {
        IdxT bus1_index = 0;
        if (branchdata.buses.contains(BranchBuses::bus1))
        {
          bus1_index = branchdata.buses.at(BranchBuses::bus1);
        }

        IdxT bus2_index = 0;
        if (branchdata.buses.contains(BranchBuses::bus2))
        {
          bus2_index = branchdata.buses.at(BranchBuses::bus2);
        }

        auto* branch = new Branch<ScalarT, IdxT>(
            getBus(bus1_index), getBus(bus2_index), branchdata);
        addComponent(branch);
      }

      // Add loads
      /// @todo Add loads to JSON parser
      for (const auto& loaddata : data.loadz)
      {
        IdxT bus_index = 0;
        if (loaddata.buses.contains(LoadZBuses::bus))
        {
          bus_index = loaddata.buses.at(LoadZBuses::bus);
        }
        auto* load = new LoadZ<ScalarT, IdxT>(getBus(bus_index), loaddata);
        addComponent(load);
      }

      // Add zip loads
      /// @todo Add zip loads to JSON parser
      for (const auto& loadzipdata : data.loadzip)
      {
        IdxT bus_index = 0;
        if (loadzipdata.buses.contains(LoadZIPBuses::bus))
        {
          bus_index = loadzipdata.buses.at(LoadZIPBuses::bus);
        }
        auto* loadzip = new LoadZIP<ScalarT, IdxT>(getBus(bus_index), loadzipdata);
        addComponent(loadzip);
      }

      // Add GENROU generators
      for (const auto& gendata : data.genrou)
      {
        IdxT bus_index = 0;
        if (gendata.buses.contains(GenrouBuses::bus))
        {
          bus_index = gendata.buses.at(GenrouBuses::bus);
        }

        auto* gen = new Genrou<ScalarT, IdxT>(getBus(bus_index), gendata);

        /// @todo Genrou (and likely other components) would need to name multiple
        /// signal inlets and outlets. For now we have only speed out and mechanical
        /// power in.
        if (gendata.signal_outputs.contains(GenrouSignalOutputs::speed))
        {
          IdxT           speed = gendata.signal_outputs.at(GenrouSignalOutputs::speed);
          constexpr auto OMEGA = GenrouInternalVariables::OMEGA;
          gen->getSignals().template assignSignalNode<OMEGA>(getSignal(speed));
        }

        if (gendata.signal_inputs.contains(GenrouSignalInputs::pmech))
        {
          IdxT           pmech = gendata.signal_inputs.at(GenrouSignalInputs::pmech);
          constexpr auto PM    = GenrouExternalVariables::PM;
          gen->getSignals().template attachSignalNode<PM>(getSignal(pmech));
        }

        if (gendata.signal_inputs.contains(GenrouSignalInputs::efd))
        {
          IdxT           efd = gendata.signal_inputs.at(GenrouSignalInputs::efd);
          constexpr auto EFD = GenrouExternalVariables::EFD;
          gen->getSignals().template attachSignalNode<EFD>(getSignal(efd));
        }

        addComponent(gen);
      }

      // Add GENSAL generators
      for (const auto& gendata : data.gensal)
      {
        IdxT bus_index = 0;
        if (gendata.buses.contains(GensalBuses::bus))
        {
          bus_index = gendata.buses.at(GensalBuses::bus);
        }

        auto* gen = new Gensal<ScalarT, IdxT>(getBus(bus_index), gendata);

        if (gendata.signal_outputs.contains(GensalSignalOutputs::speed))
        {
          IdxT           speed = gendata.signal_outputs.at(GensalSignalOutputs::speed);
          constexpr auto OMEGA = GensalInternalVariables::OMEGA;
          gen->getSignals().template assignSignalNode<OMEGA>(getSignal(speed));
        }

        if (gendata.signal_inputs.contains(GensalSignalInputs::pmech))
        {
          IdxT           pmech = gendata.signal_inputs.at(GensalSignalInputs::pmech);
          constexpr auto PM    = GensalExternalVariables::PM;
          gen->getSignals().template attachSignalNode<PM>(getSignal(pmech));
        }

        if (gendata.signal_inputs.contains(GensalSignalInputs::efd))
        {
          IdxT           efd = gendata.signal_inputs.at(GensalSignalInputs::efd);
          constexpr auto EFD = GensalExternalVariables::EFD;
          gen->getSignals().template attachSignalNode<EFD>(getSignal(efd));
        }

        addComponent(gen);
      }

      // Add classical generators
      for (const auto& gendata : data.genclassical)
      {
        IdxT bus_index = 0;
        if (gendata.buses.contains(GenClassicalBuses::bus))
        {
          bus_index = gendata.buses.at(GenClassicalBuses::bus);
        }
        auto* gen = new GenClassical<ScalarT, IdxT>(getBus(bus_index), gendata);
        addComponent(gen);
      }

      // Add REECB after its current-command and feedback producers because
      // components initialize in insertion order.
      for (const auto& reecbdata : data.reecb)
      {
        BusT* bus = nullptr;
        if (reecbdata.buses.contains(ReecbBuses::bus))
        {
          bus = getBus(reecbdata.buses.at(ReecbBuses::bus));
        }

        auto* reecb = new Reecb<ScalarT, IdxT>(bus, reecbdata);

        if (reecbdata.signal_inputs.contains(ReecbSignalInputs::pe))
        {
          const IdxT     pe = reecbdata.signal_inputs.at(ReecbSignalInputs::pe);
          constexpr auto PE = ReecbExternalVariables::PE;
          reecb->getSignals().template attachSignalNode<PE>(getSignal(pe));
        }
        if (reecbdata.signal_inputs.contains(ReecbSignalInputs::qgen))
        {
          const IdxT     qgen = reecbdata.signal_inputs.at(ReecbSignalInputs::qgen);
          constexpr auto QGEN = ReecbExternalVariables::QGEN;
          reecb->getSignals().template attachSignalNode<QGEN>(getSignal(qgen));
        }
        if (reecbdata.signal_inputs.contains(ReecbSignalInputs::qext))
        {
          const IdxT     qext = reecbdata.signal_inputs.at(ReecbSignalInputs::qext);
          constexpr auto QEXT = ReecbExternalVariables::QEXT;
          reecb->getSignals().template attachSignalNode<QEXT>(getSignal(qext));
        }
        if (reecbdata.signal_inputs.contains(ReecbSignalInputs::pfaref))
        {
          const IdxT     pfaref = reecbdata.signal_inputs.at(ReecbSignalInputs::pfaref);
          constexpr auto PFAREF = ReecbExternalVariables::PFAREF;
          reecb->getSignals().template attachSignalNode<PFAREF>(getSignal(pfaref));
        }
        if (reecbdata.signal_inputs.contains(ReecbSignalInputs::pref))
        {
          const IdxT     pref = reecbdata.signal_inputs.at(ReecbSignalInputs::pref);
          constexpr auto PREF = ReecbExternalVariables::PREF;
          reecb->getSignals().template attachSignalNode<PREF>(getSignal(pref));
        }
        if (reecbdata.signal_outputs.contains(ReecbSignalOutputs::iqcmd))
        {
          const IdxT     iqcmd = reecbdata.signal_outputs.at(ReecbSignalOutputs::iqcmd);
          constexpr auto IQCMD = ReecbInternalVariables::IQCMD;
          reecb->getSignals().template assignSignalNode<IQCMD>(getSignal(iqcmd));
        }
        if (reecbdata.signal_outputs.contains(ReecbSignalOutputs::ipcmd))
        {
          const IdxT     ipcmd = reecbdata.signal_outputs.at(ReecbSignalOutputs::ipcmd);
          constexpr auto IPCMD = ReecbInternalVariables::IPCMD;
          reecb->getSignals().template assignSignalNode<IPCMD>(getSignal(ipcmd));
        }

        addComponent(reecb);
      }

      // Add Tgov1 governors
      for (const auto& govdata : data.gov)
      {
        auto* gov = new Tgov1<ScalarT, IdxT>(govdata);

        if (govdata.signal_inputs.contains(Tgov1SignalInputs::speed))
        {
          IdxT           speed      = govdata.signal_inputs.at(Tgov1SignalInputs::speed);
          constexpr auto DELTAOMEGA = Tgov1ExternalVariables::DELTAOMEGA;
          gov->getSignals().template attachSignalNode<DELTAOMEGA>(getSignal(speed));
        }

        if (govdata.signal_outputs.contains(Tgov1SignalOutputs::pmech))
        {
          IdxT           pmech = govdata.signal_outputs.at(Tgov1SignalOutputs::pmech);
          constexpr auto PM    = Tgov1InternalVariables::PM;
          gov->getSignals().template assignSignalNode<PM>(getSignal(pmech));
        }

        addComponent(gov);
      }

      // Add GASTPTI governors
      for (const auto& gastptidata : data.gastpti)
      {
        auto* gastpti = new GastPti<ScalarT, IdxT>(gastptidata);

        if (gastptidata.signal_inputs.contains(GastPtiSignalInputs::speed))
        {
          const IdxT     speed = gastptidata.signal_inputs.at(GastPtiSignalInputs::speed);
          constexpr auto OMEGA = GastPtiExternalVariables::OMEGA;
          gastpti->getSignals().template attachSignalNode<OMEGA>(getSignal(speed));
        }

        if (gastptidata.signal_inputs.contains(GastPtiSignalInputs::pref))
        {
          const IdxT     pref = gastptidata.signal_inputs.at(GastPtiSignalInputs::pref);
          constexpr auto PREF = GastPtiExternalVariables::PREF;
          gastpti->getSignals().template attachSignalNode<PREF>(getSignal(pref));
        }

        if (gastptidata.signal_outputs.contains(GastPtiSignalOutputs::pmech))
        {
          const IdxT     pmech = gastptidata.signal_outputs.at(GastPtiSignalOutputs::pmech);
          constexpr auto PMECH = GastPtiInternalVariables::PMECH;
          gastpti->getSignals().template assignSignalNode<PMECH>(getSignal(pmech));
        }

        addComponent(gastpti);
      }

      // Add HYGOV governors
      for (const auto& hygovdata : data.hygov)
      {
        auto* hygov = new Hygov<ScalarT, IdxT>(hygovdata);

        if (hygovdata.signal_inputs.contains(HygovSignalInputs::speed))
        {
          IdxT           speed = hygovdata.signal_inputs.at(HygovSignalInputs::speed);
          constexpr auto OMEGA = HygovExternalVariables::OMEGA;
          hygov->getSignals().template attachSignalNode<OMEGA>(getSignal(speed));
        }

        if (hygovdata.signal_inputs.contains(HygovSignalInputs::pref))
        {
          IdxT           pref = hygovdata.signal_inputs.at(HygovSignalInputs::pref);
          constexpr auto PREF = HygovExternalVariables::PREF;
          hygov->getSignals().template attachSignalNode<PREF>(getSignal(pref));
        }

        if (hygovdata.signal_inputs.contains(HygovSignalInputs::paux))
        {
          IdxT           paux = hygovdata.signal_inputs.at(HygovSignalInputs::paux);
          constexpr auto PAUX = HygovExternalVariables::PAUX;
          hygov->getSignals().template attachSignalNode<PAUX>(getSignal(paux));
        }

        if (hygovdata.signal_outputs.contains(HygovSignalOutputs::pmech))
        {
          IdxT           pmech = hygovdata.signal_outputs.at(HygovSignalOutputs::pmech);
          constexpr auto PMECH = HygovInternalVariables::PMECH;
          hygov->getSignals().template assignSignalNode<PMECH>(getSignal(pmech));
        }

        addComponent(hygov);
      }

      for (const auto& excitedata : data.exciter)
      {
        IdxT bus_index = 0;
        if (excitedata.buses.contains(Ieeet1Buses::bus))
        {
          bus_index = excitedata.buses.at(Ieeet1Buses::bus);
        }

        auto* exciter = new Ieeet1<ScalarT, IdxT>(getBus(bus_index), excitedata);

        if (excitedata.signal_inputs.contains(Ieeet1SignalInputs::speed))
        {
          IdxT           speed = excitedata.signal_inputs.at(Ieeet1SignalInputs::speed);
          constexpr auto OMEGA = Ieeet1ExternalVariables::OMEGA;
          exciter->getSignals().template attachSignalNode<OMEGA>(getSignal(speed));
        }

        if (excitedata.signal_outputs.contains(Ieeet1SignalOutputs::efd))
        {
          IdxT           efd = excitedata.signal_outputs.at(Ieeet1SignalOutputs::efd);
          constexpr auto EFD = Ieeet1InternalVariables::EFD;
          exciter->getSignals().template assignSignalNode<EFD>(getSignal(efd));
        }

        if (excitedata.signal_inputs.contains(Ieeet1SignalInputs::vs))
        {
          IdxT           vs = excitedata.signal_inputs.at(Ieeet1SignalInputs::vs);
          constexpr auto VS = Ieeet1ExternalVariables::VS;
          exciter->getSignals().template attachSignalNode<VS>(getSignal(vs));
        }

        addComponent(exciter);
      }

      for (const auto& excitedata : data.esdc1a)
      {
        BusT* bus = nullptr;
        if (excitedata.buses.contains(Esdc1aBuses::bus))
        {
          bus = getBus(excitedata.buses.at(Esdc1aBuses::bus));
        }

        auto* exciter = new Esdc1a<ScalarT, IdxT>(bus, excitedata);

        if (excitedata.signal_inputs.contains(Esdc1aSignalInputs::speed))
        {
          IdxT           speed = excitedata.signal_inputs.at(Esdc1aSignalInputs::speed);
          constexpr auto OMEGA = Esdc1aExternalVariables::OMEGA;
          exciter->getSignals().template attachSignalNode<OMEGA>(getSignal(speed));
        }

        if (excitedata.signal_inputs.contains(Esdc1aSignalInputs::vref))
        {
          IdxT           vref = excitedata.signal_inputs.at(Esdc1aSignalInputs::vref);
          constexpr auto VREF = Esdc1aExternalVariables::VREF;
          exciter->getSignals().template attachSignalNode<VREF>(getSignal(vref));
        }

        if (excitedata.signal_outputs.contains(Esdc1aSignalOutputs::efd))
        {
          IdxT           efd = excitedata.signal_outputs.at(Esdc1aSignalOutputs::efd);
          constexpr auto EFD = Esdc1aInternalVariables::EFD;
          exciter->getSignals().template assignSignalNode<EFD>(getSignal(efd));
        }

        if (excitedata.signal_inputs.contains(Esdc1aSignalInputs::vs))
        {
          IdxT           vs = excitedata.signal_inputs.at(Esdc1aSignalInputs::vs);
          constexpr auto VS = Esdc1aExternalVariables::VS;
          exciter->getSignals().template attachSignalNode<VS>(getSignal(vs));
        }

        if (excitedata.signal_inputs.contains(Esdc1aSignalInputs::vuel))
        {
          IdxT           vuel = excitedata.signal_inputs.at(Esdc1aSignalInputs::vuel);
          constexpr auto VUEL = Esdc1aExternalVariables::VUEL;
          exciter->getSignals().template attachSignalNode<VUEL>(getSignal(vuel));
        }

        addComponent(exciter);
      }

      for (const auto& excitedata : data.sexspti)
      {
        IdxT bus_index = 0;
        if (excitedata.buses.contains(SexsPtiBuses::bus))
        {
          bus_index = excitedata.buses.at(SexsPtiBuses::bus);
        }

        auto* exciter = new SexsPti<ScalarT, IdxT>(getBus(bus_index), excitedata);

        if (excitedata.signal_outputs.contains(SexsPtiSignalOutputs::efd))
        {
          IdxT           efd = excitedata.signal_outputs.at(SexsPtiSignalOutputs::efd);
          constexpr auto EFD = SexsPtiInternalVariables::EFD;
          exciter->getSignals().template assignSignalNode<EFD>(getSignal(efd));
        }

        if (excitedata.signal_inputs.contains(SexsPtiSignalInputs::vs))
        {
          IdxT           vs = excitedata.signal_inputs.at(SexsPtiSignalInputs::vs);
          constexpr auto VS = SexsPtiExternalVariables::VS;
          exciter->getSignals().template attachSignalNode<VS>(getSignal(vs));
        }

        addComponent(exciter);
      }

      // Add IEEEST stabilizers
      for (const auto& stabdata : data.stabilizer)
      {
        auto* stabilizer = new Ieeest<ScalarT, IdxT>(stabdata);

        if (stabdata.signal_inputs.contains(IeeestSignalInputs::input))
        {
          IdxT           input = stabdata.signal_inputs.at(IeeestSignalInputs::input);
          constexpr auto U     = IeeestExternalVariables::U;
          stabilizer->getSignals().template attachSignalNode<U>(getSignal(input));
        }

        if (stabdata.signal_outputs.contains(IeeestSignalOutputs::output))
        {
          IdxT           output = stabdata.signal_outputs.at(IeeestSignalOutputs::output);
          constexpr auto VSS    = IeeestInternalVariables::VSS;
          stabilizer->getSignals().template assignSignalNode<VSS>(getSignal(output));
        }

        addComponent(stabilizer);
      }

      // Add REPCA plant controllers after the signal producers they read at
      // initialization
      for (const auto& repcadata : data.repca)
      {
        BusT* bus = nullptr;
        if (repcadata.buses.contains(RepcaBuses::bus))
        {
          bus = getBus(repcadata.buses.at(RepcaBuses::bus));
        }

        auto* repca = new Repca<ScalarT, IdxT>(bus, repcadata);

        if (repcadata.signal_inputs.contains(RepcaSignalInputs::ir))
        {
          const IdxT     ir = repcadata.signal_inputs.at(RepcaSignalInputs::ir);
          constexpr auto IR = RepcaExternalVariables::IR;
          repca->getSignals().template attachSignalNode<IR>(getSignal(ir));
        }
        if (repcadata.signal_inputs.contains(RepcaSignalInputs::ii))
        {
          const IdxT     ii = repcadata.signal_inputs.at(RepcaSignalInputs::ii);
          constexpr auto II = RepcaExternalVariables::II;
          repca->getSignals().template attachSignalNode<II>(getSignal(ii));
        }
        if (repcadata.signal_inputs.contains(RepcaSignalInputs::p))
        {
          const IdxT     p = repcadata.signal_inputs.at(RepcaSignalInputs::p);
          constexpr auto P = RepcaExternalVariables::P;
          repca->getSignals().template attachSignalNode<P>(getSignal(p));
        }
        if (repcadata.signal_inputs.contains(RepcaSignalInputs::q))
        {
          const IdxT     q = repcadata.signal_inputs.at(RepcaSignalInputs::q);
          constexpr auto Q = RepcaExternalVariables::Q;
          repca->getSignals().template attachSignalNode<Q>(getSignal(q));
        }
        if (repcadata.signal_inputs.contains(RepcaSignalInputs::freq))
        {
          const IdxT     freq = repcadata.signal_inputs.at(RepcaSignalInputs::freq);
          constexpr auto FREQ = RepcaExternalVariables::FREQ;
          repca->getSignals().template attachSignalNode<FREQ>(getSignal(freq));
        }
        if (repcadata.signal_inputs.contains(RepcaSignalInputs::vref))
        {
          const IdxT     vref = repcadata.signal_inputs.at(RepcaSignalInputs::vref);
          constexpr auto VREF = RepcaExternalVariables::VREF;
          repca->getSignals().template attachSignalNode<VREF>(getSignal(vref));
        }
        if (repcadata.signal_inputs.contains(RepcaSignalInputs::pref))
        {
          const IdxT     pref = repcadata.signal_inputs.at(RepcaSignalInputs::pref);
          constexpr auto PREF = RepcaExternalVariables::PREF;
          repca->getSignals().template attachSignalNode<PREF>(getSignal(pref));
        }
        if (repcadata.signal_inputs.contains(RepcaSignalInputs::qref))
        {
          const IdxT     qref = repcadata.signal_inputs.at(RepcaSignalInputs::qref);
          constexpr auto QREF = RepcaExternalVariables::QREF;
          repca->getSignals().template attachSignalNode<QREF>(getSignal(qref));
        }
        if (repcadata.signal_inputs.contains(RepcaSignalInputs::freqref))
        {
          const IdxT     freqref = repcadata.signal_inputs.at(RepcaSignalInputs::freqref);
          constexpr auto FREQREF = RepcaExternalVariables::FREQREF;
          repca->getSignals().template attachSignalNode<FREQREF>(getSignal(freqref));
        }
        if (repcadata.signal_outputs.contains(RepcaSignalOutputs::qext))
        {
          const IdxT     qext = repcadata.signal_outputs.at(RepcaSignalOutputs::qext);
          constexpr auto QEXT = RepcaInternalVariables::QEXT;
          repca->getSignals().template assignSignalNode<QEXT>(getSignal(qext));
        }
        if (repcadata.signal_outputs.contains(RepcaSignalOutputs::pext))
        {
          const IdxT     pext = repcadata.signal_outputs.at(RepcaSignalOutputs::pext);
          constexpr auto PEXT = RepcaInternalVariables::PEXT;
          repca->getSignals().template assignSignalNode<PEXT>(getSignal(pext));
        }

        addComponent(repca);
      }

      // Add constant signal sources
      for (const auto& srcdata : data.constant_source)
      {
        auto* source = new ConstantSignalSource<ScalarT, IdxT>(srcdata);

        using SignalOutputs = ConstantSignalSourceSignalOutputs;
        if (srcdata.signal_outputs.contains(SignalOutputs::sr))
        {
          IdxT           sr    = srcdata.signal_outputs.at(SignalOutputs::sr);
          constexpr auto SREAL = ConstantSignalSourceInternalVariables::SREAL;
          source->getSignals().template assignSignalNode<SREAL>(getSignal(sr));
        }
        if (srcdata.signal_outputs.contains(SignalOutputs::si))
        {
          IdxT           si    = srcdata.signal_outputs.at(SignalOutputs::si);
          constexpr auto SIMAG = ConstantSignalSourceInternalVariables::SIMAG;
          source->getSignals().template assignSignalNode<SIMAG>(getSignal(si));
        }

        addComponent(source);
      }

      // Add faults
      for (const auto& faultdata : data.bus_fault)
      {
        IdxT bus_index = 0;
        if (faultdata.buses.contains(BusFaultBuses::bus))
        {
          bus_index = faultdata.buses.at(BusFaultBuses::bus);
        }
        auto* fault = new BusFault<ScalarT, IdxT>(getBus(bus_index), faultdata);
        addFault(fault);
      }

      for (const auto& sink : data.monitor_sink)
      {
        monitor_->addSink(sink);
      }
    }

    /**
     * @brief Destructor for the system model
     *
     * If the SystemModel owns the components, it needs to delete them upon
     * destructor call.
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::~SystemModel()
    {
      if (owns_components_)
      {
        for (auto component : components_)
        {
          delete component;
        }

        for (auto bus : buses_)
        {
          delete bus;
        }

        for (auto signal : signals_)
        {
          delete signal;
        }
      }
    }

    /**
     * @brief Set component ID
     *
     * @note Should default to 0. Nested system models are not currently
     * supported.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief Allocate system storage and bind buses and components to it.
     *
     * First sum the bus and component sizes, then allocate the system vectors
     * and bind each bus and component to its portion of those vectors.
     *
     * @pre Buses and components with nonzero size are unallocated or already
     * bound to external system storage.
     *
     * @note System model composition is flat; nested systems are not supported.
     *
     * @throws std::runtime_error if storage allocation, child binding, model
     * verification, or initialization for sparse Jacobian discovery fails.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocate()
    {
      size_ = 0;

      for (const auto& bus : buses_)
      {
        size_ += bus->size();
      }

      for (const auto& component : components_)
      {
        size_ += component->size();
      }

      // Allocate global vectors
      if (!allocated_)
      {
        // Topology changes invalidate the Jacobian sparsity pattern and COO-to-CSR map.
        delete csr_jac_;
        csr_jac_ = nullptr;

        delete[] map_to_csr_;
        map_to_csr_ = nullptr;

        nnz_ = 0;
        this->allocateVectors(size_);
      }

      if (y_.getSize() != size_
          || yp_.getSize() != size_
          || f_.getSize() != size_
          || abs_tol_.getSize() != size_)
      {
        Log::error() << "SystemModel vector sizes do not match the system size\n";
        throw std::runtime_error("SystemModel allocation failed");
      }

      tag_.resize(size_);
      variable_indices_.resize(size_);
      residual_indices_.resize(size_);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      IdxT offset = 0;

      for (const auto& bus : buses_)
      {
        const int bind_status = bus->bind(y_, yp_, f_, abs_tol_, offset);
        if (bind_status != 0)
        {
          Log::error() << "Failed to bind bus vectors to system storage\n";
          throw std::runtime_error("SystemModel allocation failed");
        }

        if (bus->allocate() != 0)
        {
          Log::error() << "Failed to allocate bus data\n";
          throw std::runtime_error("SystemModel allocation failed");
        }

        for (IdxT j = 0; j < bus->size(); ++j)
        {
          bus->setVariableIndex(j, offset + j);
          bus->setResidualIndex(j, offset + j);
        }

        offset += bus->size();
      }

      for (const auto& component : components_)
      {
        const int bind_status = component->bind(y_, yp_, f_, abs_tol_, offset);
        if (bind_status != 0)
        {
          Log::error() << "Failed to bind component vectors to system storage\n";
          throw std::runtime_error("SystemModel allocation failed");
        }

        if (component->allocate() != 0)
        {
          Log::error() << "Failed to allocate component data\n";
          throw std::runtime_error("SystemModel allocation failed");
        }

        for (IdxT j = 0; j < component->size(); ++j)
        {
          component->setVariableIndex(j, offset + j);
          component->setResidualIndex(j, offset + j);
        }

        offset += component->size();
      }

      if (offset != size_)
      {
        Log::error() << "Bound vector sizes do not match the system size\n";
        throw std::runtime_error("SystemModel allocation failed");
      }

      // Verify component configuration
      int errorCount = this->verify();
      if (errorCount > 0)
      {
        Log::error() << "Component errors: " << errorCount << std::endl;
        throw std::runtime_error("SystemModel allocation failed");
      }

      // Sparse-pattern discovery requires an initialized operating point. A failed
      // initialization aborts allocation before residual/Jacobian evaluation or
      // monitor startup.
      // @todo Replace with a sparsity analysis that sets the NNZ and allocates
      // the Jacobian without needing the Jacobian values.
      if (hasJacobian())
      {
        const int status = initialize();
        if (status != 0)
        {
          Log::error() << "System model initialization failed with status "
                       << status << '\n';
          throw std::runtime_error("SystemModel allocation failed");
        }
        evaluateResidual();
        evaluateJacobian();
      }

      initializeMonitor();
      startMonitor();

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Verify all components are configured correctly
     *
     * This method accumulates and returns the number of errors given by
     * components. It should return 0 when all is well.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::verify() const
    {
      int ret = 0;

      // Verify components
      for (const auto& component : components_)
      {
        ret += component->verify();
      }

      return ret;
    }

    /**
     * @brief Initialize buses first, then all the other components.
     *
     * @pre All buses and components must be allocated at this point.
     * @pre Bus variables are written before component variables in the
     * system variable vector.
     *
     * Buses must be initialized before other components, because other
     * components may write to buses during the initialization.
     *
     * Also, generators may write to control devices (e.g. governors,
     * exciters, etc.) during the initialization.
     *
     * @todo Implement writing to system vectors in a thread-safe way.
     *
     * @note Currently assuming each component stores variables contiguously in memory and
     * that these are simply concateneted in the global system.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::initialize()
    {
      int status = 0;

      for (const auto& bus : buses_)
      {
        status += bus->initialize();
      }

      for (const auto& component : components_)
      {
        status += component->initialize();
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return status;
    }

    /**
     * @brief Add monitors from buses and components and start monitor
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::initializeMonitor()
    {
      for (const auto* bus : buses_)
      {
        auto* mon = bus->getMonitor();
        if (mon && !mon->empty())
        {
          monitor_->addMonitor(mon);
        }
      }

      for (const auto* component : components_)
      {
        auto* mon = component->getMonitor();
        if (mon && !mon->empty())
        {
          monitor_->addMonitor(mon);
        }
      }
    }

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::startMonitor()
    {
      monitor_->start();
    }

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::stopMonitor()
    {
      monitor_->stop();
    }

    template <typename scalar_type, typename index_type>
    bool SystemModel<scalar_type, index_type>::monitoring() const
    {
      return !monitor_->empty();
    }

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::printMonitoredVariables() const
    {
      monitor_->print();
    }

    /**
     * @todo Tagging differential variables
     *
     * Identify what variables in the system of differential-algebraic
     * equations are differential variables, i.e. their derivatives
     * appear in the equations.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::tagDifferentiable()
    {
      // Set initial values for global solution vectors
      for (const auto& bus : buses_)
      {
        bus->tagDifferentiable();
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          tag_[bus->getVariableIndex(j)] = bus->tag()[j];
        }
      }

      for (const auto& component : components_)
      {
        component->tagDifferentiable();
        for (IdxT j = 0; j < component->size(); ++j)
        {
          tag_[component->getVariableIndex(j)] = component->tag()[j];
        }
      }

      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      for (const auto& bus : buses_)
      {
        bus->setAbsoluteTolerance(rel_tol);
      }

      for (const auto& component : components_)
      {
        component->setAbsoluteTolerance(rel_tol);
      }

      abs_tol_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Compute system residual vector
     *
     * Buses and components read and write their bound system-vector slices
     * directly.
     *
     * @warning Residuals must be computed for buses, before component
     * residuals are computed. Buses own residuals for currents
     * Ir and Ii, but the contributions to these residuals come
     * from components. Buses assign their residual values, while components
     * add to those values by in-place adition. This is why (for now) bus
     * residuals need to be computed first.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateResidual()
    {
      for (const auto& bus : buses_)
      {
        bus->evaluateResidual();
      }

      for (const auto& component : components_)
      {
        component->evaluateResidual();
      }

      f_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Evaluate system Jacobian.
     *
     * First, initialize bus Jacobians to 0.
     * Then, evaluate component Jacobians (internal block and bus Jacobian contributions).
     * Once component Jacobians are evaluated, store the result in the system Jacobian.
     * Finally, store bus Jacobians into the system Jacobian after all component have added their
     * contributions.
     *
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      // Initialize bus Jacobians
      for (const auto& bus : buses_)
      {
        bus->evaluateJacobian();
      }

      // Evaluate component Jacobians, including contribution to the bus Jacobians
      for (const auto& component : components_)
      {
        component->evaluateJacobian();
      }

      // Build or update system CSR Jacobian
      if (csr_jac_ == nullptr)
      {
        // Count the number of non-zeros
        IdxT nnz_dup = 0;
        for (const auto& component : components_)
        {
          auto component_jacobian = component->getCooJacobian();

          if (component_jacobian != nullptr)
          {
            nnz_dup += component_jacobian->getNnz();
          }
          else
          {
            Log::warning() << "A component has returned a nullptr Jacobian.\n";
          }
        }

        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getCooJacobian();

          if (bus_jacobian != nullptr)
          {
            nnz_dup += bus_jacobian->getNnz();
          }
          else
          {
            Log::warning() << "A bus has returned a nullptr Jacobian.\n";
          }
        }

        // Allocate COO triplet arrays (we own these until we hand off to CsrMatrix)
        IdxT*  rows_dup = new IdxT[static_cast<size_t>(nnz_dup)];
        IdxT*  cols_dup = new IdxT[static_cast<size_t>(nnz_dup)];
        RealT* vals_dup = new RealT[static_cast<size_t>(nnz_dup)];

        IdxT counter = 0;
        for (const auto& component : components_)
        {
          auto component_jacobian = component->getCooJacobian();

          if (component_jacobian != nullptr)
          {
            const IdxT*  rows    = component_jacobian->getRowData();
            const IdxT*  columns = component_jacobian->getColData();
            const RealT* values  = component_jacobian->getValues();
            for (IdxT i = 0; i < component_jacobian->getNnz(); ++i)
            {
              rows_dup[counter] = rows[i];
              cols_dup[counter] = columns[i];
              vals_dup[counter] = values[i];
              counter++;
            }
          }
          else
          {
            Log::warning() << "A component has returned a nullptr Jacobian.\n";
          }
        }

        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getCooJacobian();

          if (bus_jacobian != nullptr)
          {
            const IdxT*  rows    = bus_jacobian->getRowData();
            const IdxT*  columns = bus_jacobian->getColData();
            const RealT* values  = bus_jacobian->getValues();
            for (IdxT i = 0; i < bus_jacobian->getNnz(); ++i)
            {
              rows_dup[counter] = rows[i];
              cols_dup[counter] = columns[i];
              vals_dup[counter] = values[i];
              counter++;
            }
          }
          else
          {
            Log::warning() << "A bus has returned a nullptr Jacobian.\n";
          }
        }

        // Build the system COO Jacobian
        CooMatrixT jac(size_, size_, nnz_dup, &rows_dup, &cols_dup, &vals_dup);

        // Populate CSR data with sort and deduplicate
        IdxT* row_ptrs = jac.getCsrRowData();

        // Deduplicated nnz
        nnz_ = jac.getNnz();

        // Allocate cols/vals with deduplicated nnz
        IdxT*  cols = new IdxT[static_cast<size_t>(nnz_)];
        RealT* vals = new RealT[static_cast<size_t>(nnz_)];

        std::copy(jac.getColData(), jac.getColData() + nnz_, cols);
        std::copy(jac.getValues(), jac.getValues() + nnz_, vals);

        // Create the CSR Jacobian
        csr_jac_ = new CsrMatrixT(size_, size_, nnz_, &row_ptrs, &cols, &vals);

        const IdxT* map_to_sorted = jac.getMapToSorted();
        const IdxT* map_to_dedup  = jac.getMapToDeduplicated();

        // Build a mappping from original COO index to CSR index
        map_to_csr_ = new IdxT[static_cast<size_t>(nnz_dup)];
        for (IdxT i = 0; i < nnz_dup; ++i)
        {
          map_to_csr_[map_to_sorted[i]] = map_to_dedup[i];
        }
      }
      else
      {
        // Zero out values
        RealT* vals = csr_jac_->getValues();
        for (IdxT i = 0; i < csr_jac_->getNnz(); ++i)
        {
          vals[i] = 0.0;
        }

        // Update CSR values from component and bus Jacobians
        IdxT counter = 0;
        for (const auto& component : components_)
        {
          auto component_jacobian = component->getCooJacobian();

          if (component_jacobian != nullptr)
          {
            const RealT* values = component_jacobian->getValues();
            for (IdxT i = 0; i < component_jacobian->getNnz(); ++i)
            {
              vals[map_to_csr_[counter]] += values[i];
              counter++;
            }
          }
        }

        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getCooJacobian();

          if (bus_jacobian != nullptr)
          {
            const RealT* values = bus_jacobian->getValues();
            for (IdxT i = 0; i < bus_jacobian->getNnz(); ++i)
            {
              vals[map_to_csr_[counter]] += values[i];
              counter++;
            }
          }
        }
      }

      // std::cout << "System Jacobian\n";
      // csr_jac_->print(std::cout);

      return 0;
    }

    /**
     * @brief Update time
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::updateTime(RealT t, RealT a)
    {
      time_  = t;
      alpha_ = a;
      for (const auto& component : components_)
      {
        component->updateTime(t, a);
      }
    }

    /**
     * @brief Add bus
     *
     * Add bus at the end of the bus array and map bus ID with GridKit's ID for the bus
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addBus(BusT* bus)
    {
      IdxT gridkit_bus_id                = static_cast<IdxT>(buses_.size());
      gridkit_bus_indices_[bus->busID()] = gridkit_bus_id;
      buses_.push_back(bus);
      allocated_ = false;
    }

    /**
     * @brief Add signal
     *
     * Add signal at the end of the signals array and map signal ID with
     * GridKit's ID for the signal
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addSignal(SignalT* signal)
    {
      IdxT gridkit_signal_id                      = static_cast<IdxT>(signals_.size());
      gridkit_signal_indices_[signal->signalId()] = gridkit_signal_id;
      signals_.push_back(signal);
    }

    /**
     * @brief Add component
     *
     * Add component at the end of the components array and set GridKit's component ID
     *
     * @todo: No integer user-facing component_id for now, but we could map GridKit's
     * component ID to the disambiguation_string
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addComponent(ComponentT* component)
    {
      IdxT gridkit_component_id = static_cast<IdxT>(components_.size());
      component->setGridKitComponentID(gridkit_component_id);
      component->setSystemBase(this->freq_system_base_,
                               this->va_system_base_);
      components_.push_back(component);
      allocated_ = false;
    }

    /**
     * @brief Add fault
     *
     * The fault is added to the components array, and we keep a map to its
     * location, so it can easily be accessed.
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addFault(ComponentT* component)
    {
      IdxT gridkit_component_id                = static_cast<IdxT>(components_.size());
      IdxT gridkit_fault_id                    = static_cast<IdxT>(gridkit_fault_indices_.size());
      gridkit_fault_indices_[gridkit_fault_id] = gridkit_component_id;
      addComponent(component);
    }

    /**
     * @brief Set system bases and propagate them to existing components.
     *
     * @param[in] freq_system_base - System frequency base in Hz.
     * @param[in] va_system_base - System power base in VA.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::setSystemBase(
        RealT freq_system_base, RealT va_system_base)
    {
      ComponentT::setSystemBase(freq_system_base, va_system_base);

      for (auto* component : components_)
      {
        component->setSystemBase(freq_system_base, va_system_base);
      }
    }

    /**
     * @brief Return pointer to a bus
     *
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::BusT*
    SystemModel<scalar_type, index_type>::getBus(IdxT bus_id)
    {
      // Should fail if user-provided bus_id is incorrect
      IdxT gridkit_bus_id = gridkit_bus_indices_.at(bus_id);
      assert((buses_[gridkit_bus_id])->busID() == bus_id);
      return buses_[gridkit_bus_id];
    }

    /**
     * @brief Return pointer to a signal
     *
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SignalT*
    SystemModel<scalar_type, index_type>::getSignal(IdxT signal_id)
    {
      // Should fail if user-provided signal_id is incorrect
      IdxT gridkit_signal_id = gridkit_signal_indices_.at(signal_id);
      assert((signals_[gridkit_signal_id])->signalId() == signal_id);
      return signals_[gridkit_signal_id];
    }

    /**
     * @brief Return pointer to a component
     *
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::ComponentT*
    SystemModel<scalar_type, index_type>::getComponent(IdxT gridkit_component_id)
    {
      // gridkit_component_id_ is set by System model and guarantied to be unique
      return components_[gridkit_component_id];
    }

    /**
     * @brief Return pointer to a bus fault model
     *
     * This function is used to provide easier access to setting and
     * clearing faults from the SystemModel interface.
     *
     */
    template <typename scalar_type, typename index_type>
    BusFault<scalar_type, index_type>*
    SystemModel<scalar_type, index_type>::getBusFault(IdxT fault_id)
    {
      IdxT component_id = gridkit_fault_indices_.at(fault_id);
      return dynamic_cast<BusFault<ScalarT, IdxT>*>(components_[component_id]);
    }

  } // namespace PhasorDynamics
} // namespace GridKit
