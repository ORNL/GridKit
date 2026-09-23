/**
 * @file StateData.hpp
 * @brief Operating state shared by all model families. See `STATE.md`.
 */

#pragma once

#include <cstddef>
#include <filesystem>
#include <istream>
#include <map>
#include <optional>
#include <ostream>
#include <string>
#include <utility>

namespace GridKit
{
  namespace Model
  {
    /// Entries of one bus or device
    struct StateRecord
    {
      std::map<std::string, double> values; ///< Numeric entries such as `vr` or `ir`
      std::map<std::string, bool>   flags;  ///< Boolean entries such as `online` or `open`

      double value(const std::string& key, double fallback) const;
      bool   flag(const std::string& key, bool fallback) const;
    };

    /// Metadata of a state file
    struct StateHeader
    {
      std::optional<unsigned int> version;
      std::optional<double>       time;
      std::optional<std::string>  created;
      std::optional<std::string>  description;
    };

    /**
     * @brief Operating state of a case
     *
     * Currents are positive into the bus and use the system base.
     */
    struct StateData
    {
      StateHeader                        header;
      std::map<std::string, StateRecord> buses;   ///< Keyed by `busKey(number)`
      std::map<std::string, StateRecord> devices; ///< Keyed by the case `id`

      const StateRecord* bus(std::size_t number) const;
      const StateRecord* device(const std::string& id) const;
    };

    /// Key of bus `number` in `StateData::buses`
    std::string busKey(std::size_t number);

    /**
     * @brief Keys of the real and imaginary current at `terminal` of a device
     * with `terminal_count` terminals
     */
    std::pair<std::string, std::string> currentKeys(std::size_t terminal, std::size_t terminal_count);

    /**
     * @brief Power into bus `number` at `terminal` of device `id`
     *
     * @return false if the state has no current for the terminal or no
     * voltage for the bus
     */
    bool terminalPower(const StateData&   state,
                       const std::string& id,
                       std::size_t        number,
                       std::size_t        terminal,
                       std::size_t        terminal_count,
                       double&            p,
                       double&            q);

    /**
     * @brief Set the current of device `id` at `terminal` that injects `p`
     * and `q` into bus `number`
     *
     * @pre The state has a nonzero voltage for bus `number`.
     */
    void setTerminalCurrent(StateData&         state,
                            const std::string& id,
                            std::size_t        number,
                            std::size_t        terminal,
                            std::size_t        terminal_count,
                            double             p,
                            double             q);

    StateData parseStateData(std::istream& stream);
    StateData parseStateData(const std::filesystem::path& file);
    void      writeStateData(const StateData& state, std::ostream& stream);
    void      writeStateData(const StateData& state, const std::filesystem::path& file);
  } // namespace Model
} // namespace GridKit
