

#ifndef _CIRCCOMP_HPP_
#define _CIRCCOMP_HPP_

#include <ModelEvaluatorImpl.hpp>
#include <PowerSystemData.hpp>
#include <map>
#include <set>


namespace ModelLib
{
    /*!
     * @brief Declaration of a passive CircuitComponent class.
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

        inline std::vector<ScalarT> parkTransformMatrix(ScalarT angle)
        {
            std::vector<ScalarT> result(9);
            ScalarT anpim = angle - (2.0/3.0)*M_PI;
            ScalarT anpip = angle + (2.0/3.0)*M_PI;
            result[0] = (2.0/3.0)*cos(angle);
            result[1] = (2.0/3.0)*cos(anpim);
            result[2] = (2.0/3.0)*cos(anpip);
            result[3] = (2.0/3.0)*sin(angle);
            result[4] = (2.0/3.0)*sin(anpim);
            result[5] = (2.0/3.0)*sin(anpip);
            result[6] = 1.0/3.0;
            result[7] = 1.0/3.0;
            result[8] = 1.0/3.0;
            return result;
        }

        inline std::vector<ScalarT> parkTransformMatrixDerivative(ScalarT angle)
        {
            std::vector<ScalarT> result(9);
            ScalarT anpim = angle - (2.0/3.0)*M_PI;
            ScalarT anpip = angle + (2.0/3.0)*M_PI;
            result[0] = (2.0/3.0)*sin(angle);
            result[1] = (2.0/3.0)*sin(anpim);
            result[2] = (2.0/3.0)*sin(anpip);
            result[3] = (2.0/3.0)*-cos(angle);
            result[4] = (2.0/3.0)*-cos(anpim);
            result[5] = (2.0/3.0)*-cos(anpip);
            result[6] = 0;
            result[7] = 0;
            result[8] = 0;
            return result;
        }

        inline std::vector<ScalarT> inverseParkTransformMatrix(ScalarT angle)
        {
            std::vector<ScalarT> result(9);
            ScalarT anpim = angle - (2.0/3.0)*M_PI;
            ScalarT anpip = angle + (2.0/3.0)*M_PI;
            result[0] = cos(angle);
            result[1] = sin(angle);
            result[2] = 1.0;
            result[3] = cos(anpim);
            result[4] = sin(anpim);
            result[5] = 1.0;
            result[6] = cos(anpip);
            result[7] = sin(anpip);
            result[8] = 1.0;
            return result;
        }

        inline std::vector<ScalarT> inverseParkTransformMatrixDerivative(ScalarT angle)
        {
            std::vector<ScalarT> result(9);
            ScalarT anpim = angle - (2.0/3.0)*M_PI;
            ScalarT anpip = angle + (2.0/3.0)*M_PI;
            result[0] = sin(angle);
            result[1] = -cos(angle);
            result[2] = 0.0;
            result[3] = sin(anpim);
            result[4] = -cos(anpim);
            result[5] = 0.0;
            result[6] = sin(anpip);
            result[7] = -cos(anpip);
            result[8] = 0.0;
            return result;
        }

    protected:
        size_t n_extern_;
        size_t n_intern_;
        std::set<IdxT> extern_indices_;
        //@todo may want to replace the mapping of connection_nodes to Node objects instead of IdxT. Allows for container free setup
        std::map<IdxT, IdxT> connection_nodes_;

    };

    
}

#endif
