/**
 * @file matpower_parser.hpp
 *
 * Standalone MATPOWER v2 .m file parser for use in pf-solver/.
 * Kept separate from GridKit source so it can be modified freely without
 * touching the GridKit tree. Reads only the fields needed for PF:
 *   mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost
 * Unknown mpc.* fields (mpc.bus_name, mpc.gentype, mpc.genfuel, etc.) are
 * silently skipped. Extra columns beyond the expected count are ignored.
 *
 * Output: populates GridKit::PowerFlowData::SystemModelData<RealT, IdxT>.
 */
#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <GridKit/Model/PowerFlow/PowerFlowData.hpp>

namespace UqPfSolver
{
  using namespace GridKit::PowerFlowData;

  // Known mpc.* matrix fields that we parse.
  static const std::string KNOWN_FIELDS[] = {"bus", "gen", "branch", "gencost"};

  static inline void ltrim(std::string& s)
  {
    auto it = s.begin();
    while (it != s.end() && std::isspace((unsigned char)*it)) ++it;
    s.erase(s.begin(), it);
  }

  static inline void rtrim(std::string& s)
  {
    while (!s.empty() && std::isspace((unsigned char)s.back())) s.pop_back();
  }

  // Strip everything from '%' to end of line.
  static inline void strip_comments(std::string& s)
  {
    auto pos = s.find('%');
    if (pos != std::string::npos) s.erase(pos);
  }

  // Extract mpc.<field> name from a line like "mpc.bus = [" or "mpc.baseMVA = 100;".
  // Returns empty string if line does not match.
  static inline std::string get_mpc_field(const std::string& line)
  {
    auto dot = line.find("mpc.");
    if (dot == std::string::npos) return "";
    auto start = dot + 4;
    auto end   = start;
    while (end < line.size() &&
           (std::isalnum((unsigned char)line[end]) || line[end] == '_'))
      ++end;
    return line.substr(start, end - start);
  }

  // Skip lines until the line containing "];" is consumed.
  static inline void skip_matrix(std::istream& is)
  {
    std::string line;
    while (std::getline(is, line))
    {
      strip_comments(line);
      if (line.find("];") != std::string::npos) return;
    }
  }

  // Parse baseMVA from "mpc.baseMVA = 100;" already trimmed.
  template <typename RealT, typename IdxT>
  static void parse_basemva(const std::string& line,
                             SystemModelData<RealT, IdxT>& mp)
  {
    auto eq = line.find('=');
    if (eq == std::string::npos) return;
    std::string val = line.substr(eq + 1);
    strip_comments(val);
    ltrim(val); rtrim(val);
    // Remove trailing semicolon
    if (!val.empty() && val.back() == ';') val.pop_back();
    mp.baseMVA = std::stod(val);
  }

  // Read bus matrix rows. MATPOWER bus columns (1-indexed):
  //  1=bus_i 2=type 3=Pd 4=Qd 5=Gs 6=Bs 7=area 8=Vm 9=Va 10=baseKV 11=zone 12=Vmax 13=Vmin
  // Extra columns are ignored.
  template <typename RealT, typename IdxT>
  static void parse_bus_matrix(std::istream& is,
                                SystemModelData<RealT, IdxT>& mp)
  {
    std::string line;
    while (std::getline(is, line))
    {
      strip_comments(line);
      ltrim(line); rtrim(line);
      if (line.empty()) continue;
      if (line.find("];") != std::string::npos) return;

      std::istringstream ss(line);
      BusData<RealT, IdxT>  bd{};
      LoadData<RealT, IdxT> ld{};
      RealT dummy;
      // cols: bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
      if (!(ss >> bd.bus_i >> bd.type >> ld.Pd >> ld.Qd
               >> bd.Gs >> bd.Bs >> bd.area >> bd.Vm >> bd.Va
               >> bd.baseKV >> bd.zone >> bd.Vmax >> bd.Vmin))
        continue; // skip malformed or empty rows
      // MATPOWER Va is in degrees; GridKit Bus residuals use radians internally
      bd.Va *= M_PI / 180.0;
      ld.bus_i = bd.bus_i;
      // ignore any extra columns
      mp.bus.push_back(std::move(bd));
      mp.load.push_back(std::move(ld));
    }
  }

  // Read gen matrix rows. MATPOWER gen columns (1-indexed):
  //  1=bus 2=Pg 3=Qg 4=Qmax 5=Qmin 6=Vg 7=mBase 8=status 9=Pmax 10=Pmin
  //  11=Pc1 12=Pc2 13=Qc1min 14=Qc1max 15=Qc2min 16=Qc2max 17=ramp_agc
  //  18=ramp_10 19=ramp_30 20=ramp_q 21=apf
  // Extended MATPOWER cases (e.g. ACTIVSg) add extra columns -- ignored.
  template <typename RealT, typename IdxT>
  static void parse_gen_matrix(std::istream& is,
                                SystemModelData<RealT, IdxT>& mp)
  {
    std::string line;
    while (std::getline(is, line))
    {
      strip_comments(line);
      ltrim(line); rtrim(line);
      if (line.empty()) continue;
      if (line.find("];") != std::string::npos) return;

      std::istringstream ss(line);
      GenData<RealT, IdxT> gd{};
      if (!(ss >> gd.bus >> gd.Pg >> gd.Qg >> gd.Qmax >> gd.Qmin
               >> gd.Vg >> gd.mBase >> gd.status >> gd.Pmax >> gd.Pmin
               >> gd.Pc1 >> gd.Pc2 >> gd.Qc1min >> gd.Qc1max
               >> gd.Qc2min >> gd.Qc2max >> gd.ramp_agc >> gd.ramp_10
               >> gd.ramp_30 >> gd.ramp_q >> gd.apf))
        continue;
      // ignore any extra columns
      mp.gen.push_back(gd);
    }
  }

  // Read branch matrix rows. MATPOWER branch columns (1-indexed):
  //  1=fbus 2=tbus 3=r 4=x 5=b 6=rateA 7=rateB 8=rateC 9=ratio 10=angle
  //  11=status 12=angmin 13=angmax
  // Extended cases add extra columns -- ignored.
  template <typename RealT, typename IdxT>
  static void parse_branch_matrix(std::istream& is,
                                   SystemModelData<RealT, IdxT>& mp)
  {
    std::string line;
    while (std::getline(is, line))
    {
      strip_comments(line);
      ltrim(line); rtrim(line);
      if (line.empty()) continue;
      if (line.find("];") != std::string::npos) return;

      std::istringstream ss(line);
      BranchData<RealT, IdxT> bd{};
      if (!(ss >> bd.fbus >> bd.tbus >> bd.r >> bd.x >> bd.b
               >> bd.rateA >> bd.rateB >> bd.rateC >> bd.ratio
               >> bd.angle >> bd.status >> bd.angmin >> bd.angmax))
        continue;
      // ignore any extra columns
      mp.branch.push_back(bd);
    }
  }

  // Read gencost matrix rows (not used by KINSOL PF solve but parsed for
  // completeness so the loop below can skip it cleanly).
  template <typename RealT, typename IdxT>
  static void parse_gencost_matrix(std::istream& is,
                                    SystemModelData<RealT, IdxT>& mp)
  {
    std::string line;
    while (std::getline(is, line))
    {
      strip_comments(line);
      ltrim(line); rtrim(line);
      if (line.empty()) continue;
      if (line.find("];") != std::string::npos) return;

      std::istringstream ss(line);
      GenCostData<RealT, IdxT> gc{};
      if (!(ss >> gc.kind >> gc.startup >> gc.shutdown >> gc.n))
        continue;
      RealT v;
      while (ss >> v) gc.rest.push_back(v);
      mp.gencost.push_back(gc);
    }
  }

  /**
   * Parse a MATPOWER v2 .m file into SystemModelData.
   * Handles:
   *   - fields with underscores (mpc.bus_name, mpc.gen_type, etc.) -- skipped
   *   - extra columns in gen/branch rows (ACTIVSg extended format) -- ignored
   *   - cell arrays ({ }) -- skipped
   */
  template <typename RealT = double, typename IdxT = size_t>
  void readMatPowerFile(SystemModelData<RealT, IdxT>& mp,
                        const std::string& filename)
  {
    std::ifstream ifs{filename};
    if (!ifs)
      throw std::runtime_error("cannot open: " + filename);

    std::string line;
    while (std::getline(ifs, line))
    {
      strip_comments(line);
      ltrim(line); rtrim(line);
      if (line.empty()) continue;
      if (line.find("function") != std::string::npos) continue;
      if (line.find("mpc.") == std::string::npos) continue;

      std::string field = get_mpc_field(line);
      if (field.empty()) continue;

      if (field == "baseMVA")
      {
        parse_basemva(line, mp);
      }
      else if (field == "version")
      {
        // parse version string -- not critical, skip
      }
      else if (field == "bus")
      {
        // line is "mpc.bus = [" -- matrix starts on next line
        parse_bus_matrix(ifs, mp);
      }
      else if (field == "gen")
      {
        parse_gen_matrix(ifs, mp);
      }
      else if (field == "branch")
      {
        parse_branch_matrix(ifs, mp);
      }
      else if (field == "gencost")
      {
        parse_gencost_matrix(ifs, mp);
      }
      else
      {
        // Unknown field (bus_name, gentype, genfuel, ...).
        // If it opens a matrix "[" or cell array "{", consume until "];" or "};"
        if (line.find('[') != std::string::npos ||
            line.find('{') != std::string::npos)
        {
          skip_matrix(ifs);
        }
        // scalar fields: nothing to consume beyond this line
      }
    }

    // Normalize MW/MVAr quantities to per-unit.
    // GridKit's PF model (Branch, Generator, Load) all work in pu;
    // MATPOWER .m files store Pd/Qd/Pg/Qg in MW/MVAr.
    const RealT inv_base = (mp.baseMVA > 0.0) ? 1.0 / mp.baseMVA : 1.0;
    for (auto& ld : mp.load)
    {
      ld.Pd *= inv_base;
      ld.Qd *= inv_base;
    }
    for (auto& gd : mp.gen)
    {
      gd.Pg  *= inv_base;
      gd.Qg  *= inv_base;
      gd.Qmax *= inv_base;
      gd.Qmin *= inv_base;
      gd.Pmax *= inv_base;
      gd.Pmin *= inv_base;
    }

    std::cerr << "Parsed (local parser): baseMVA=" << mp.baseMVA
              << "  buses=" << mp.bus.size()
              << "  gens=" << mp.gen.size()
              << "  branches=" << mp.branch.size() << "\n";
  }

  /**
   * @brief Write a solved .m file by patching Vm/Va in mpc.bus of the source file.
   *
   * Reads the original .m text line by line. When inside the mpc.bus matrix,
   * replaces columns 8 (Vm) and 9 (Va, degrees) in each data row with the
   * converged PF solution stored in @p bus_solution (keyed by bus_i).
   * Everything else (mpc.gen, mpc.branch, comments, other fields) is copied
   * through unchanged.
   *
   * @param src_path   Path to the original (input) .m file.
   * @param dst_path   Path to write the solved (output) .m file.
   * @param bus_solution  Map from bus_i -> (Vm_pu, Va_deg) from the PF solution.
   */
  inline void write_solved_m(
      const std::string& src_path,
      const std::string& dst_path,
      const std::map<size_t, std::pair<double, double>>& bus_solution)
  {
    std::ifstream ifs{src_path};
    if (!ifs) throw std::runtime_error("write_solved_m: cannot open input: " + src_path);

    std::ofstream ofs{dst_path};
    if (!ofs) throw std::runtime_error("write_solved_m: cannot open output: " + dst_path);

    bool in_bus_matrix = false;
    std::string line;
    while (std::getline(ifs, line))
    {
      if (!in_bus_matrix)
      {
        // Detect start of mpc.bus matrix.
        std::string stripped = line;
        strip_comments(stripped);
        ltrim(stripped); rtrim(stripped);
        std::string field = get_mpc_field(stripped);
        if (field == "bus" && stripped.find('[') != std::string::npos)
        {
          in_bus_matrix = true;
          ofs << line << "\n";
          continue;
        }
        ofs << line << "\n";
      }
      else
      {
        // Inside mpc.bus: patch or pass through each row.
        std::string stripped = line;
        strip_comments(stripped);
        ltrim(stripped); rtrim(stripped);

        // End of matrix.
        if (stripped.find("];") != std::string::npos ||
            stripped == "]")
        {
          in_bus_matrix = false;
          ofs << line << "\n";
          continue;
        }

        // Skip blank or comment-only lines.
        if (stripped.empty())
        {
          ofs << line << "\n";
          continue;
        }

        // Parse enough columns to get bus_i (col 0), then locate and replace
        // cols 7 (Vm) and 8 (Va) by token position in the original line.
        std::istringstream ss(stripped);
        size_t bus_i = 0;
        if (!(ss >> bus_i))
        {
          ofs << line << "\n";
          continue;
        }

        auto it = bus_solution.find(bus_i);
        if (it == bus_solution.end())
        {
          // Bus not in solution (e.g. isolated bus) — pass through unchanged.
          ofs << line << "\n";
          continue;
        }

        double vm_new = it->second.first;
        double va_new = it->second.second;  // degrees

        // Tokenize the original line to find and replace cols 7 and 8.
        // Preserve the original whitespace structure by replacing tokens in-place.
        // We work on the stripped (comment-removed) content but emit the original
        // prefix whitespace.

        // Find leading whitespace prefix from original line.
        std::string prefix;
        for (char c : line)
        {
          if (std::isspace((unsigned char)c)) prefix += c;
          else break;
        }

        // Split stripped into tokens.
        std::vector<std::string> tokens;
        {
          std::istringstream ts(stripped);
          std::string tok;
          while (ts >> tok) tokens.push_back(tok);
        }

        // Replace token[7] (Vm) and token[8] (Va).
        if (tokens.size() > 8)
        {
          char buf[64];
          std::snprintf(buf, sizeof(buf), "%.8f", vm_new);
          tokens[7] = buf;
          std::snprintf(buf, sizeof(buf), "%.8f", va_new);
          tokens[8] = buf;
        }

        // Reassemble with single spaces and original indentation.
        ofs << prefix;
        for (size_t i = 0; i < tokens.size(); ++i)
        {
          if (i > 0) ofs << "\t";
          ofs << tokens[i];
        }
        // Preserve trailing semicolon if the original had one.
        if (line.find(';') != std::string::npos &&
            (tokens.empty() || tokens.back() != ";"))
          ofs << ";";
        ofs << "\n";
      }
    }
    std::cerr << "Wrote solved .m: " << dst_path << "\n";
  }

} // namespace UqPfSolver
