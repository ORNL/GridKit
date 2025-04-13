/**
 * @file Genrou.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a GENROU generator model.
 *
 *
 */

 #include "ClassicalGen.hpp"

 #include <cmath>
 #include <iostream>
 
 #include <Model/PhasorDynamics/Bus/Bus.hpp>
 
 #define _USE_MATH_DEFINES
 
 namespace GridKit
 {
   namespace PhasorDynamics
   {
     /*!
      * @brief Constructor for a pi-model branch
      *
      * Arguments passed to ModelEvaluatorImpl:
      * - Number of equations = 0
      * - Number of independent variables = 0
      * - Number of quadratures = 0
      * - Number of optimization parameters = 0
      */
     template <class ScalarT, typename IdxT>
     ClassicalGen<ScalarT, IdxT>::ClassicalGen(bus_type* bus, int unit_id)
       : bus_(bus),
         busID_(0),
         unit_id_(unit_id),
         H_(3.),
         D_(0.),
         Ra_(0.),
         Xdp_(.5),
         pmech_(0.2),
         ep_(0.2)
     {
       size_ = 5;
       setDerivedParams();
 
       // Temporary, to eliminate compiler warnings
       (void) busID_;
       (void) unit_id_;
     }
 
     /*!
      * @brief Constructor for a pi-model branch
      *
      * Arguments passed to ModelEvaluatorImpl:
      * - Number of equations = 0
      * - Number of independent variables = 0
      * - Number of quadratures = 0
      * - Number of optimization parameters = 0
      */
     template <class ScalarT, typename IdxT>
     ClassicalGen<ScalarT, IdxT>::ClassicalGen(bus_type* bus,
                                   int       unit_id,
                                   real_type H,
                                   real_type D,
                                   real_type Ra,
                                   real_type Xdp,
                                   real_type pmech,
                                   real_type ep)
                                   
       : bus_(bus),
         busID_(0),
         unit_id_(unit_id),
         H_(H),
         D_(D),
         Ra_(Ra),
         Xdp_(Xdp),
         pmech_(pmech),
         ep_(ep)
     {
       size_ = 5;
       setDerivedParams();
     }
 
     // /**
     //  * @brief Destroy the Genrou
     //  *
     //  * @tparam ScalarT
     //  * @tparam IdxT
     //  */
     // template <class ScalarT, typename IdxT>
     // Genrou<ScalarT, IdxT>::~Genrou()
     // {
     //   // std::cout << "Destroy Genrou..." << std::endl;
     // }
 
     /*!
      * @brief allocate method computes sparsity pattern of the Jacobian.
      */
     template <class ScalarT, typename IdxT>
     int ClassicalGen<ScalarT, IdxT>::allocate()
     {
       f_.resize(size_);
       y_.resize(size_);
       yp_.resize(size_);
       tag_.resize(size_);
       fB_.resize(size_);
       yB_.resize(size_);
       ypB_.resize(size_);
       return 0;
     }
 
     /**
      * Initialization of the branch model
      *
      */
     template <class ScalarT, typename IdxT>
     int ClassicalGen<ScalarT, IdxT>::initialize()
     {
       
       return 0;
     }
 
     /**
      * \brief Identify differential variables.
      */
     template <class ScalarT, typename IdxT>
     int ClassicalGen<ScalarT, IdxT>::tagDifferentiable()
     {

       return 0;
     }
 
     /**
      * \brief Residual contribution of the branch is pushed to the
      * two terminal buses.
      *
      */
     template <class ScalarT, typename IdxT>
     int ClassicalGen<ScalarT, IdxT>::evaluateResidual()
     {
       /* Read variables */
       ScalarT delta  = y_[0];
       ScalarT omega  = y_[1];
       ScalarT telec  = y_[2];
       ScalarT ir     = y_[3];
       ScalarT ii     = y_[4];

       /* Read derivatives */
       ScalarT delta_dot = yp_[0];
       ScalarT omega_dot = yp_[1];

       /* 6 ClassicalGen differential equations */
       f_[0] = delta_dot - omega * (2 * M_PI * 60);
       f_[1] = omega_dot - (1.0 / (2 * H_)) * ((pmech_ - D_ * omega) / (1 + omega) - telec);
       
       /* 11 ClassicalGen algebraic equations */
       f_[2] = telec - (1.0/(1.0 + omega))*(g*ep_*ep_ - ep_*(cos(delta)*(g*Vr() - b*Vi()) + sin(delta)*(b*Vr() + g*Vi())));

       f_[3] = ir + g*Vr() - b * Vi()  - ep_*(g*cos(delta) -b*sin(delta));
       f_[4] = ii + b*Vr() +  g * Vi() - ep_*(b*cos(delta) + g*sin(delta));

       Ir() += - (g*Vr() - b * Vi()  - ep_*(g*cos(delta) -b*sin(delta)));
       Ii() += - (b*Vr() +  g * Vi() - ep_*(b*cos(delta) + g*sin(delta)));

       return 0;
     }
 
     /**
      * @brief Jacobian evaluation not implemented yet
      *
      * @tparam ScalarT - scalar data type
      * @tparam IdxT    - matrix index data type
      * @return int - error code, 0 = success
      */
     template <class ScalarT, typename IdxT>
     int ClassicalGen<ScalarT, IdxT>::evaluateJacobian()
     {
       return 0;
     }
 
     /**
      * @brief Integrand (objective) evaluation not implemented yet
      *
      * @tparam ScalarT - scalar data type
      * @tparam IdxT    - matrix index data type
      * @return int - error code, 0 = success
      */
     template <class ScalarT, typename IdxT>
     int ClassicalGen<ScalarT, IdxT>::evaluateIntegrand()
     {
       // std::cout << "Evaluate Integrand for Genrou..." << std::endl;
       return 0;
     }
 
     /**
      * @brief Adjoint initialization not implemented yet
      *
      * @tparam ScalarT - scalar data type
      * @tparam IdxT    - matrix index data type
      * @return int - error code, 0 = success
      */
     template <class ScalarT, typename IdxT>
     int ClassicalGen<ScalarT, IdxT>::initializeAdjoint()
     {
       // std::cout << "Initialize adjoint for Genrou..." << std::endl;
       return 0;
     }
 
     /**
      * @brief Adjoint residual evaluation not implemented yet
      *
      * @tparam ScalarT - scalar data type
      * @tparam IdxT    - matrix index data type
      * @return int - error code, 0 = success
      */
     template <class ScalarT, typename IdxT>
     int ClassicalGen<ScalarT, IdxT>::evaluateAdjointResidual()
     {
       // std::cout << "Evaluate adjoint residual for Genrou..." << std::endl;
       return 0;
     }
 
     /**
      * @brief Adjoint integrand (objective) evaluation not implemented yet
      *
      * @tparam ScalarT - scalar data type
      * @tparam IdxT    - matrix index data type
      * @return int - error code, 0 = success
      */
     template <class ScalarT, typename IdxT>
     int ClassicalGen<ScalarT, IdxT>::evaluateAdjointIntegrand()
     {
       // std::cout << "Evaluate adjoint Integrand for Genrou..." << std::endl;
       return 0;
     }
 
     template <class ScalarT, typename IdxT>
     void ClassicalGen<ScalarT, IdxT>::setDerivedParams()
     {
       g   = Ra_ / (Ra_ * Ra_ + Xdp_ * Xdp_);
       b   = Xdp_ / (Ra_ * Ra_ + Xdp_ * Xdp_);
     }
 
     // Available template instantiations
     template class ClassicalGen<double, long int>;
     template class ClassicalGen<double, size_t>;
 
   } // namespace PhasorDynamics
 } // namespace GridKit
 