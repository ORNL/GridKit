/**
 * @file SignalNodeSet.hpp
 * @author superwhiskers <whiskerdev@protonmail.com>
 * @brief Collection of signal nodes.
 */

#pragma once

#include <cassert>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Manage a collection of signal nodes specified for a system model
     *
     * This class manages the memory and access of signal nodes for a
     * SystemModel. A SignalNode can be accessed either by id or name.
     */
    template <typename scalar_type, typename index_type>
    class SignalNodeSet
    {
    public:
      using ScalarT         = scalar_type;
      using IdxT            = index_type;
      using SignalNodeT     = SignalNode<ScalarT, IdxT>;
      using RealT           = typename SignalNodeT::RealT;
      using SignalNodeDataT = SignalNodeData<RealT, IdxT>;

      SignalNodeSet()                                = default;
      SignalNodeSet(const SignalNodeSet&)            = delete;
      SignalNodeSet& operator=(const SignalNodeSet&) = delete;

      ~SignalNodeSet()
      {
        for (auto node : signal_nodes_)
        {
          delete node;
        }
      }

      /**
       * @brief Create a new node based on the input data
       */
      void add(const SignalNodeDataT& node_data)
      {
        auto node = new SignalNodeT(node_data);
        signal_nodes_.push_back(node);
        signal_node_id_map_[node_data.signal_id] = node;
        signal_node_name_map_[node_data.name]    = node;
      }

      /**
       * @brief Access a SignalNode using its id
       *
       * @throws std::runtime_error if no node has this ID
       */
      SignalNodeT* operator[](IdxT signal_id)
      {
        if (!signal_node_id_map_.contains(signal_id))
        {
          Log::error() << "Unknown signal ID " << signal_id << "\n";
          throw std::runtime_error("SignalNodeSet lookup failed");
        }

        auto* node = signal_node_id_map_.at(signal_id);
        assert(node->signalId() == signal_id);
        return node;
      }

      /**
       * @brief Return the SignalNode with this id, or nullptr if there is none
       */
      SignalNodeT* find(IdxT signal_id) const
      {
        SignalNodeT* node = nullptr;
        const auto   it   = signal_node_id_map_.find(signal_id);
        if (it != signal_node_id_map_.end())
        {
          node = it->second;
        }
        return node;
      }

      /**
       * @brief Access a SignalNode using its name
       */
      SignalNodeT* operator[](const std::string& signal_name)
      {
        return signal_node_name_map_.at(signal_name);
      }

      /**
       * @brief Signal nodes in the order they were added
       */
      const std::vector<SignalNodeT*>& nodes() const
      {
        return signal_nodes_;
      }

    private:
      /// Vector of nodes
      std::vector<SignalNodeT*>                     signal_nodes_;
      /// Map of ids to nodes
      std::unordered_map<IdxT, SignalNodeT*>        signal_node_id_map_;
      /// Map of names to nodes
      std::unordered_map<std::string, SignalNodeT*> signal_node_name_map_;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
