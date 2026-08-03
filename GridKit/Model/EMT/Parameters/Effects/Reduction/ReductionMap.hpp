/**
 * @file ReductionMap.hpp
 *
 * @brief Conductor-to-terminal incidence for bundle merging and
 * grounded-conductor elimination.
 *
 */

#pragma once

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/ConductorData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      /**
       * @brief Conductor-to-terminal incidence derived from the parsed
       *        phase and circuit labels.
       *
       * Conductors sharing a phase and circuit form one electrical
       * terminal; conductors labeled with the ground phase belong to no
       * terminal. Terminals are ordered by first appearance among the
       * conductors, so the mapping is deterministic for a given line
       * description.
       */
      template <typename scalar_type, typename index_type>
      struct ReductionMap
      {
        using ScalarT = scalar_type;
        using IdxT    = index_type;

        /// Terminal index of a grounded conductor.
        static constexpr IdxT kGrounded = INVALID_INDEX<IdxT>;

        IdxT conductors{0};
        IdxT terminals{0};

        /// terminal[c] for bundle members; kGrounded for grounded ones.
        std::vector<IdxT> terminal;

        /// members[t]: conductor indices of terminal t, conductor order.
        std::vector<std::vector<IdxT>> members;

        /// Every conductor is its own terminal and none is grounded, so
        /// the reduction would be the identity.
        bool identity() const
        {
          return terminals == conductors;
        }

        static ReductionMap
        fromConductors(const ConductorData<ScalarT, IdxT>& data)
        {
          ReductionMap map;
          map.conductors = data.K;
          map.terminal.assign(static_cast<size_t>(data.K), kGrounded);

          std::vector<std::pair<std::string, IdxT>> keys;
          for (IdxT c = 0; c < data.K; ++c)
          {
            const auto& phase = data.phase[static_cast<size_t>(c)];
            if (phase == "g")
            {
              continue;
            }
            const IdxT circuit = static_cast<size_t>(c) < data.circuit.size()
                                     ? data.circuit[static_cast<size_t>(c)]
                                     : IdxT{1};

            IdxT index = static_cast<IdxT>(keys.size());
            for (size_t k = 0; k < keys.size(); ++k)
            {
              if (keys[k].first == phase && keys[k].second == circuit)
              {
                index = static_cast<IdxT>(k);
                break;
              }
            }
            if (index == static_cast<IdxT>(keys.size()))
            {
              keys.emplace_back(phase, circuit);
              map.members.emplace_back();
            }
            map.terminal[static_cast<size_t>(c)] = index;
            map.members[static_cast<size_t>(index)].push_back(c);
          }

          map.terminals = static_cast<IdxT>(keys.size());
          if (map.conductors > 0 && map.terminals == 0)
          {
            throw std::runtime_error(
                "Line reduction requires at least one ungrounded conductor");
          }
          return map;
        }
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
