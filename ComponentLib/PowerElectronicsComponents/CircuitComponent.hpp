

#ifndef _CIRCCOMP_HPP_
#define _CIRCCOMP_HPP_

#include <ModelEvaluatorImpl.hpp>
#include <PowerSystemData.hpp>
#include <map>
#include <set>


namespace ModelLib
{
    /*!
     * @brief Declaration of a CircuitComponent class.
     *
     */
    template  <class ScalarT, typename IdxT>
    class CircuitComponent : public ModelEvaluatorImpl<ScalarT, IdxT>
    {

    public:


        void updateTime(ScalarT t, ScalarT a)
        {
            this->time_ = t;
            this->alpha_ = a;
        }
        
        bool hasJacobian() { return true;}

        size_t getExternSize()
        {
            return n_extern_;
        }
        
        size_t getInternalSize()
        {
            return this->n_intern_;
        }

        std::set<size_t> getExternIndices()
        {
            return this->extern_indices_;
        }

        /**
         * @brief Create the mappings from local to global indexes
         * 
         * @param local_index 
         * @param global_index 
         * @return int
         */
        int setExternalConnectionNodes(IdxT local_index, IdxT global_index)
        {
            connection_nodes_[local_index] = global_index;
            return 0;
        }

        /**
         * @brief Given the location of value in the local vector map to global index
         * 
         * f(local_index) = global_index
         *
         * @param local_index index of local value in vector
         * @return IdxT Index of the same value in the global vector
         */
        IdxT getNodeConnection(IdxT local_index)
        {
            return connection_nodes_.at(local_index);
        }

    protected:
        size_t n_extern_;
        size_t n_intern_;
        std::set<IdxT> extern_indices_;
        ///@todo may want to replace the mapping of connection_nodes to Node objects instead of IdxT. Allows for container free setup
        std::map<IdxT, IdxT> connection_nodes_;

    };

    
}

#endif
