
#include <iostream>
#include <cmath>
#include <vector>
#include "TransmissionLine.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a constant TransmissionLine model
 *
 * Calls default ModelEvaluatorImpl constructor.
 * 
 * This is the Medium distance form with the use of the admittance matrix.
 * Since the line is of medium length then there is no real part for shunt admittance
 */

template <class ScalarT, typename IdxT>
TransmissionLine<ScalarT, IdxT>::TransmissionLine(IdxT id, ScalarT R,ScalarT X, ScalarT B)
  : R_(R), X_(X), B_(B) 
{
    // internals [Iret1, Iimt1, Iret2, Iimt2]
    // externals [Vre11, Vim11, Vre12, Vim12, Vre21, Vim21, Vre22, Vim22]
	this->size_ = 12;
    this->n_intern = 4;
    this->n_extern = 8;
    this->extern_indices = {0,1,2,3,4,5,6,7};
    this->idc_ = id;

    ScalarT magImpendence = 1 / (R_*R_ + X_*X_);
    YReMat_ = magImpendence * R_;
    YImMatOff_ = magImpendence * X_;
    YImMatDi_ = B_ / (2.0) -  YImMatOff_;
}

template <class ScalarT, typename IdxT>
TransmissionLine<ScalarT, IdxT>::~TransmissionLine()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int TransmissionLine<ScalarT, IdxT>::allocate()
{
    this->y_.resize(this->size_);
    this->yp_.resize(this->size_);
    this->f_.resize(this->size_);
    
    return 0;
}

/**
 * Initialization of the grid model
 */
template <class ScalarT, typename IdxT>
int TransmissionLine<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int TransmissionLine<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Evaluate residual of transmission line
 *
 * The complex admittance matrix is:
 * [[ Y/2 + 1/Z, -1/Z];
 *  [ -1/Z, Y/2 + 1/Z]] = 
 * [[R/|Z|, -R/|Z|];
 *  [-R/|Z|, R/|Z|]] + 
 * i [[B/2 - X/|Z|, X/|Z|];
 *    [X/|Z|, B/2 - X/|Z|]]
 * = Dre + i Dim
 * 
 * Then 
 *  Ire = Dre Vre - Dim Vim
 *  Iim = Dre Vim + Dim Vre
 * 
 * To express this for Modified Nodal Analysis the Voltages of the admittance matrix are put into voltage drops
 */
template <class ScalarT, typename IdxT>
int TransmissionLine<ScalarT, IdxT>::evaluateResidual()
{
	//input
    this->f_[0] = y_[8] ;
    this->f_[1] = y_[9] ;

    this->f_[2] = y_[10] ;
    this->f_[3] = y_[11] ;
    //ouput
    this->f_[4] = -y_[8] ;
    this->f_[5] = -y_[9] ;

    this->f_[6] = -y_[10] ;
    this->f_[7] = -y_[11] ;

    //Voltage drop accross terminals
    ScalarT V1re = y_[0] - y_[4];
    ScalarT V1im = y_[1] - y_[5];
    ScalarT V2re = y_[2] - y_[6];
    ScalarT V2im = y_[3] - y_[7];

    //Internal variables
    //row 1
    this->f_[8] = YReMat_ * (V1re - V2re) - (YImMatDi_ * V1im + YImMatOff_ * V2im) - y_[8] ;
    this->f_[9] = YReMat_ * (V1im - V2im) + (YImMatDi_ * V1re + YImMatOff_ * V2re) - y_[9]  ;

    //row2
    this->f_[10] = YReMat_ * (V2re - V1re) - (YImMatOff_ * V1im + YImMatDi_ * V2im) - y_[10];
    this->f_[11] = YReMat_ * (V2im - V1im) + (YImMatOff_ * V1re + YImMatDi_ * V2re) - y_[11];


    return 0;
}

/**
 * @brief Generate Jacobian for Transmission Line
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return int 
 */
template <class ScalarT, typename IdxT>
int TransmissionLine<ScalarT, IdxT>::evaluateJacobian()
{
    
    //Create dF/dy
    std::vector<IdxT> rtemp{0,1,2,3,4,5,6,7,8,9,10,11};
    std::vector<IdxT> ctemp{8,9,10,11,8,9,10,11,8,9,10,11};
    std::vector<ScalarT> vals{1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
    this->J_.setValues(rtemp, ctemp, vals);

    
    std::vector<IdxT> ccord{0,1,2,3,4,5,6,7};

    std::vector<IdxT> rcord(ccord.size(),8);
    vals = {YReMat_, -YImMatDi_ ,-YReMat_, -YImMatOff_,-YReMat_, YImMatDi_ ,YReMat_, YImMatOff_};
    this->J_.setValues(rtemp, ctemp, vals);

    
    std::fill(rcord.begin(), rcord.end(), 9);
    vals = {YImMatDi_ ,YReMat_, YImMatOff_, -YReMat_,-YImMatDi_ ,-YReMat_, -YImMatOff_, YReMat_};
    this->J_.setValues(rtemp, ctemp, vals);

    
    std::fill(rcord.begin(), rcord.end(), 10);
    vals = {-YReMat_, -YImMatDi_ ,YReMat_, -YImMatOff_,YReMat_, YImMatDi_ ,-YReMat_, YImMatOff_};
    this->J_.setValues(rtemp, ctemp, vals);

    
    std::fill(rcord.begin(), rcord.end(), 11);
    vals = {YImMatDi_ ,-YReMat_, YImMatOff_, YReMat_,-YImMatDi_ ,YReMat_, -YImMatOff_, -YReMat_};
    this->J_.setValues(rtemp, ctemp, vals);

    return 0;
}

template <class ScalarT, typename IdxT>
int TransmissionLine<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int TransmissionLine<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int TransmissionLine<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int TransmissionLine<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}





// Available template instantiations
template class TransmissionLine<double, long int>;
template class TransmissionLine<double, size_t>;


} //namespace ModelLib

