/*!
    \file VariableImplementation.hpp
*/

#pragma once

namespace GridKit
{
namespace Sparse
{
    /*!
        \brief Returns the list of derivatives.
    */
    const Variable::DependencyMap& Variable::getDependencies() const
    {
        assert(dependencies_ != 0);
        return *dependencies_;
    }
    

    /*!
        \brief Registers a variable as an unknown of the system and 
        adds a pointer to the global @a x vector.
    */
    void Variable::registerVariable(std::vector<Variable*>& x, 
                                    const size_t& offset)
    {
        setVariableNumber(offset);  // define global variable number
        setFixed(false);            // not a constant

        x[offset] = this;
    }
    

    /*!
        \brief Adds all dependencies of v to *this.
    */
    void Variable::addDependencies(const Variable& v)
    {
        for (auto& p : *(v.dependencies_))
            (*dependencies_)[p.first] = p.second;
    }
    

    /*!
        \brief Multiplies each partial derivative of @a this by @a c. 
    */
    void Variable::scaleDependencies(double c)
    {
        for (auto& p : *dependencies_)
            (*dependencies_)[p.first] *= c;
    }
    
    
    /*!
        \brief Prints the value and input set of the variable.
    */
    void Variable::print(std::ostream& os) const
    {
        os << value_;
        
        if (is_fixed_)
        {
            os << " (fixed)";
            return;
        }
        
        if (variable_number_ != INVALID_VAR_NUMBER)
        {
            os << " (variable " << variable_number_ << ")";
            return;
        }
        
        if (dependencies_!= NULL && !dependencies_->empty())
        {
            os << " dependencies: [ ";
            for (auto& p : *dependencies_)
                os << "(" << p.first << ", " << p.second << ") ";
            os << "]";
        }
    }
    
    
    //------------------------------------
    // compound assignment operators
    //------------------------------------
    
    /*!
        \brief Compound addition-assignment operator. Right hand 
        side is a built-in double type.
    */
    Variable& Variable::operator+=(const double& rhs)
    {
        value_ += rhs;
        return *this;
    }
    
    /*!
        \brief Compound addition-assignment operator. Right hand side 
        is Variable type.
    */
    Variable& Variable::operator+=(const Variable& rhs)
    {
        // compute partial derivatives of *this
        for (auto& p : *(rhs.dependencies_))
            (*dependencies_)[p.first] += (p.second);

        // compute value of *this
        value_ += rhs.value_;

        return *this;
    }
    
    // -=
    /*!
        \brief Compound subtraction-assignment operator. Right hand 
        side is a built-in double type.
    */
    Variable& Variable::operator-=(const double& rhs)
    {
        value_ -= rhs;
        return *this;
    }
    
    /*!
        \brief Compound subtraction-assignment operator. Right hand 
        side is a Variable type.
    */
    Variable& Variable::operator-=(const Variable& rhs)
    {
        // compute partial derivatives of *this
        for (auto& p : *(rhs.dependencies_))
            (*dependencies_)[p.first] -= (p.second);

        // compute value of *this
        value_ -= rhs.value_;

        return *this;
    }
    
    // *=
    /*!
        \brief Compound multiplication-assignment operator. Right hand 
        side is a built-in double type.
    */
    Variable& Variable::operator*=(const double& rhs)
    {
        // Compute derivatives of this
        scaleDependencies(rhs);

        // compute value
        value_ *= rhs;

        return *this;
    }
    
    /*!
        \brief Compound multiplication-assignment operator. Right 
        hand side is a Variable type.
    */
    Variable& Variable::operator*=(const Variable& rhs)
    {
        // derivation by parts of @ *this
        scaleDependencies(rhs.value_);

        // compute partial derivatives of rhs and add them to *this
        for (auto& p : *(rhs.dependencies_))
            (*dependencies_)[p.first] += (p.second * value_);

        // compute value of this
        value_ *= rhs.value_;

        return *this;
    }
    
    // /=
    /*!
        \brief Compound division-assignment operator. Right hand side 
        is a built-in double type.
    */
    Variable& Variable::operator/=(const double& rhs)
    {
        double inverseRhs = 1.0/rhs;

        // compute derivatives of @a *this
        scaleDependencies(inverseRhs);

        value_ *= inverseRhs;
        return *this;
    }
    
    /*!
        \brief Compound division-assignment operator. Right hand 
        side is a Variable type.
    */
    Variable& Variable::operator/=(const Variable& rhs)
    {
        double inverseRhs        = 1.0/rhs.value_;
        double inverseRhsSq = inverseRhs * inverseRhs;

        // derivation by parts of @a *this
        scaleDependencies(inverseRhs);
        for (auto& p : *(rhs.dependencies_))
            (*dependencies_)[p.first] -= (p.second * value_ * inverseRhsSq);

        // compute value of this
        value_ *= inverseRhs;

        return *this;
    }
    
} // namespace Sparse
} // namespace GridKit

