install(FILES DependencyTracking/Variable.hpp DependencyTracking/VariableImplementation.hpp
              DependencyTracking/VariableOperators.hpp
        DESTINATION include/GridKit/AutomaticDifferentiation/DependencyTracking)

if(GRIDKIT_ENABLE_ENZYME)
  install(FILES Enzyme/EnzymeDefinitions.hpp
          DESTINATION include/GridKit/AutomaticDifferentiation/Enzyme)

  install(
    FILES Enzyme/Optimization/ConstraintJacobian.hpp
          Enzyme/Optimization/LagrangianHessian.hpp
          Enzyme/Optimization/ModelWrappers.hpp
          Enzyme/Optimization/ObjectiveGradient.hpp
          Enzyme/Optimization/SparseDerivatives.hpp
    DESTINATION include/GridKit/AutomaticDifferentiation/Enzyme/Optimization)
endif()
