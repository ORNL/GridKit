

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

		size_t getExternSize()
		{
			return this->n_extern;
		}
		
		size_t getInternalSize()
		{
			return this->n_intern;
		}

		std::set<size_t> getExternIndices()
		{
			return this->extern_indices;
		}

		bool setExternalConnectionNodes(size_t index, IdxT id)
		{
			this->connection_nodes[index] = id;
			return true;
		}

		IdxT getNodeConnection(size_t index)
		{
			return this->connection_nodes.at(index);
		}

	protected:
		size_t n_extern;
		size_t n_intern;
		std::set<size_t> extern_indices;
		std::map<size_t, IdxT> connection_nodes;

    };

	
}

#endif
