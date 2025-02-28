/**
    \file Variable.hpp
*/
#pragma once

#include <map>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>

#include <cassert>
#include <limits>
#include <iostream>
#include <iomanip>
#include <string>

namespace GridKit
{
namespace Sparse
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
              is_fixed_(false),
              dependencies_(new DependencyMap)
        {
        }

        /**
            @brief Constructor which initializes the value.
        */
        explicit Variable(double value)
            : value_(value),
              variable_number_(INVALID_VAR_NUMBER),
              is_fixed_(false),
              dependencies_(new DependencyMap)
        {
        }

        /**
            @brief Constructor which initializes the value and variable 
            number.

        */
        Variable(double value, size_t variable_number)
            : value_(value),
              variable_number_(variable_number),
              is_fixed_(false),
              dependencies_(new DependencyMap)
        {
            (*dependencies_)[variable_number_] = 1.0;
        }
    
        /**
            @brief Copy constructor.
        */
        Variable(const Variable& v)
            : value_(v.value_),
              variable_number_(INVALID_VAR_NUMBER),
              is_fixed_(false),
              dependencies_(new DependencyMap(*v.dependencies_))
        {
        }
    
        /**
            @brief Destructor deletes the dependency map.
        */
        ~Variable()
        {
            delete dependencies_;
        }
    
        /**
            @brief Assignment operator. Assigning double value to 
            Variable removes its dependencies. Use only if you know 
            what you are doing.
        */
        Variable& operator=(const double& rhs)
        {
            value_ = rhs;
        
            dependencies_->clear();
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
            if (this == &rhs)  // self-assignment
                return *this;
        
            // set value from rhs
            value_ = rhs.value_;

            // if rhs is a constant so is *this
            setFixed(rhs.is_fixed());

            // set dependencies from rhs
            dependencies_->clear(); // clear map just in case
            addDependencies(rhs);   // use only dependencies from the rhs
        
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
            return (*dependencies_)[i];
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
            dependencies_->clear();
            variable_number_ = variable_number;
            (*dependencies_)[variable_number_] = 1.0;
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
                                     const size_t& offset);

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
        
        
    private:
        double value_;           ///< Value of the variable.        
        size_t variable_number_;  ///< Independent variable ID
        bool is_fixed_;           ///< Constant parameter flag.
        
        mutable DependencyMap* dependencies_;
        static const size_t INVALID_VAR_NUMBER = static_cast<const size_t>(-1);
    };
    
    //------------------------------------
    // non-member operators and functions
    //------------------------------------
    
    // unary -
    inline const Variable operator-(const Variable& v);
    
    // +
    inline const Variable operator+(const Variable& lhs, const Variable& rhs);
    inline const Variable operator+(const Variable& lhs, const double& rhs);
    inline const Variable operator+(const double& lhs, const Variable& rhs);
    
    // -
    inline const Variable operator-(const Variable& lhs, const Variable& rhs);
    inline const Variable operator-(const Variable& lhs, const double& rhs);
    inline const Variable operator-(const double& lhs, const Variable& rhs);
    
    // *
    inline const Variable operator*(const Variable& lhs, const Variable& rhs);
    inline const Variable operator*(const Variable& lhs, const double& rhs);
    inline const Variable operator*(const double& lhs, const Variable& rhs);
    
    // /
    inline const Variable operator/(const Variable& lhs, const Variable& rhs);
    inline const Variable operator/(const Variable& lhs, const double& rhs);
    inline const Variable operator/(const double& lhs, const Variable& rhs);
    
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
    inline Variable& operator<<(Variable& u, const Variable& v);
    
    inline std::istream& operator>>(std::istream& is, Variable& v);
    

} // namespace Sparse
} // namespace GridKit

#include "VariableImplementation.hpp"
#include "VariableOperators.hpp"

