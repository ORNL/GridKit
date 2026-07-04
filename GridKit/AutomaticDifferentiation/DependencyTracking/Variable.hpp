/**
    \file Variable.hpp
*/
#pragma once

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace DependencyTracking
  {

    /**
        @brief The Variable class is used to store the unknowns of the
        system and define abstract elementary algebra.

        \author Stefan Klus
        \author Slaven Peles
    */
    class Variable
    {
    public:
      /**
          @brief Default constructor.
      */
      Variable()
        : value_(0.0),
          variable_number_(INVALID_VAR_NUMBER),
          is_fixed_(false)
      {
      }

      /**
          @brief Constructor which initializes the value.
      */
      explicit Variable(double value)
        : value_(value),
          variable_number_(INVALID_VAR_NUMBER),
          is_fixed_(false)
      {
      }

      /**
          @brief Constructor which initializes the value and variable
          number.

      */
      Variable(double value, size_t variable_number)
        : value_(value),
          variable_number_(variable_number),
          is_fixed_(false)
      {
        dependencies_[variable_number_] = 1.0;
      }

      /**
          @brief Copy constructor.
      */
      Variable(const Variable& v)
        : value_(v.value_),
          variable_number_(INVALID_VAR_NUMBER),
          is_fixed_(false),
          dependencies_(v.dependencies_)
      {
      }

      /**
          @brief Move constructor.

          Unlike the copy constructor, this relocates @a v rather than
          creating a new, detached temporary from it, so (unlike copy)
          it carries over `variable_number_`/`is_fixed_` as-is instead
          of resetting them. The dependency map is relocated via
          `std::map`'s own move constructor (O(1)) instead of being
          deep-copied.
      */
      Variable(Variable&& v) noexcept
        : value_(v.value_),
          variable_number_(v.variable_number_),
          is_fixed_(v.is_fixed_),
          dependencies_(std::move(v.dependencies_))
      {
      }

      /**
          @brief Destructor.
      */
      ~Variable() = default;

      /**
          @brief Assignment operator. Assigning double value to
          Variable removes its dependencies. Use only if you know
          what you are doing.
      */
      Variable& operator=(const double& rhs)
      {
        value_ = rhs;

        dependencies_.clear();
        return *this;
      }

      /**
          @brief Assignment operator.

          This operator:
          - assigns value from the right hand side
          - leaves variable ID unchanged
          - clears any existing and adds new dependencies from rhs
      */
      Variable& operator=(const Variable& rhs)
      {
        if (this == &rhs) // self-assignment
          return *this;

        // set value from rhs
        value_ = rhs.value_;

        // if rhs is a constant so is *this
        setFixed(rhs.is_fixed());

        // set dependencies from rhs
        dependencies_.clear(); // clear map just in case
        addDependencies(rhs);  // use only dependencies from the rhs

        return *this;
      }

      /**
          @brief Move assignment operator.

          Same contract as the copy assignment operator (value and
          dependencies come from rhs, variable ID of *this is left
          unchanged), but relocates rhs's dependency map via
          `std::map`'s own move assignment instead of deep-copying its
          entries.
      */
      Variable& operator=(Variable&& rhs) noexcept
      {
        if (this == &rhs) // self-assignment
          return *this;

        // set value from rhs
        value_ = rhs.value_;

        // if rhs is a constant so is *this
        setFixed(rhs.is_fixed());

        // relocate dependencies from rhs instead of copying them
        dependencies_ = std::move(rhs.dependencies_);

        return *this;
      }

      /**
          @brief Operator () returns the value of a variable.

          This is just short notation to avoid using
          getValue and setValue.
      */
      double& operator()()
      {
        return value_;
      }

      /**
          @brief Operator() returns the value of a variable
          (const version).

          This is just short notation to avoid using getValue.
      */
      const double& operator()() const
      {
        return value_;
      }

      /**
          @brief Return the current value of the variable.
      */
      double getValue() const
      {
        return value_;
      }

      /**
          @brief Overwrite the current value of the variable.
      */
      void setValue(double value)
      {
        value_ = value;
      }

      /**
          @brief Return derivative of *this with respect to
          dependency i.
      */
      double der(size_t i) const
      {
        return dependencies_[i];
      }

      /**
          @brief Returns the variable number.

          This number is assigned to state variables (variables
          updated directly by the solver) only.
      */
      size_t getVariableNumber() const
      {
        return variable_number_;
      }

      /**
          @brief Sets the variable number.
      */
      void setVariableNumber(size_t variable_number)
      {
        dependencies_.clear();
        variable_number_                  = variable_number;
        dependencies_[variable_number_] = 1.0;
      }

      /**
          @brief Checks whether the variable was registered as
          an unknown of the system.

          INVALID_VAR_NUMBER is used to mark parameters and
          temporary variables
      */
      bool isRegistered() const
      {
        return variable_number_ != INVALID_VAR_NUMBER;
      }

      /**
          @brief Checks whether the variable is fixed or not.
      */
      bool is_fixed() const
      {
        return is_fixed_;
      }

      /**
          @brief Turns variable into parameter, or vice versa.
       */
      void setFixed(bool b = false)
      {
        is_fixed_ = b;
      }

      // get the 'input set' of a variable
      using DependencyMap = std::map<size_t, double>;
      inline const DependencyMap& getDependencies() const;

      // set as the independent state variable and assign ID to it
      inline void registerVariable(std::vector<Variable*>& x,
                                   const size_t&           offset);

      // adds all dependencies of v to *this
      inline void addDependencies(const Variable& v);

      // scale dependencies (derivatives) by scalar @a c.
      inline void scaleDependencies(double c);

      // print to output stream
      inline void print(std::ostream& os) const;

      // +=
      inline Variable& operator+=(const double& rhs);
      inline Variable& operator+=(const Variable& rhs);

      // -=
      inline Variable& operator-=(const double& rhs);
      inline Variable& operator-=(const Variable& rhs);

      // *=
      inline Variable& operator*=(const double& rhs);
      inline Variable& operator*=(const Variable& rhs);

      // /=
      inline Variable& operator/=(const double& rhs);
      inline Variable& operator/=(const Variable& rhs);

      // conversion operator
      explicit inline operator double() const
      {
        return value_;
      }

    private:
      double value_;           ///< Value of the variable.
      size_t variable_number_; ///< Independent variable ID
      bool   is_fixed_;        ///< Constant parameter flag.

      mutable DependencyMap dependencies_;
      static const size_t    INVALID_VAR_NUMBER = INVALID_INDEX<size_t>;
    };

    //------------------------------------
    // non-member operators and functions
    //------------------------------------

    // unary -
    inline Variable operator-(const Variable& v);

    // +
    inline Variable operator+(const Variable& lhs, const Variable& rhs);
    inline Variable operator+(const Variable& lhs, const double& rhs);
    inline Variable operator+(const double& lhs, const Variable& rhs);

    // -
    inline Variable operator-(const Variable& lhs, const Variable& rhs);
    inline Variable operator-(const Variable& lhs, const double& rhs);
    inline Variable operator-(const double& lhs, const Variable& rhs);

    // *
    inline Variable operator*(const Variable& lhs, const Variable& rhs);
    inline Variable operator*(const Variable& lhs, const double& rhs);
    inline Variable operator*(const double& lhs, const Variable& rhs);

    // /
    inline Variable operator/(const Variable& lhs, const Variable& rhs);
    inline Variable operator/(const Variable& lhs, const double& rhs);
    inline Variable operator/(const double& lhs, const Variable& rhs);

    // ==
    inline bool operator==(const Variable& lhs, const Variable& rhs);
    inline bool operator==(const Variable& lhs, const double& rhs);
    inline bool operator==(const double& lhs, const Variable& rhs);

    // !=
    inline bool operator!=(const Variable& lhs, const Variable& rhs);
    inline bool operator!=(const Variable& lhs, const double& rhs);
    inline bool operator!=(const double& lhs, const Variable& rhs);

    // <
    inline bool operator<(const Variable& lhs, const Variable& rhs);
    inline bool operator<(const Variable& lhs, const double& rhs);
    inline bool operator<(const double& lhs, const Variable& rhs);

    // >
    inline bool operator>(const Variable& lhs, const Variable& rhs);
    inline bool operator>(const Variable& lhs, const double& rhs);
    inline bool operator>(const double& lhs, const Variable& rhs);

    // <=
    inline bool operator<=(const Variable& lhs, const Variable& rhs);
    inline bool operator<=(const Variable& lhs, const double& rhs);
    inline bool operator<=(const double& lhs, const Variable& rhs);

    // >=
    inline bool operator>=(const Variable& lhs, const Variable& rhs);
    inline bool operator>=(const Variable& lhs, const double& rhs);
    inline bool operator>=(const double& lhs, const Variable& rhs);

    inline std::ostream& operator<<(std::ostream& os, const Variable& v);
    inline Variable&     operator<<(Variable& u, const Variable& v);

    inline std::istream& operator>>(std::istream& is, Variable& v);

  } // namespace DependencyTracking
} // namespace GridKit

namespace GridKit
{
  template <>
  class ScalarTraits<DependencyTracking::Variable>
  {
  public:
    typedef double RealT;
    typedef double NormT;
  };

} // namespace GridKit

#include "VariableImplementation.hpp"
#include "VariableOperators.hpp"
