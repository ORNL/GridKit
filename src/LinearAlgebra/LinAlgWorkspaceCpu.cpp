#include "LinAlgWorkspaceCpu.hpp"

namespace GridKit
{
  namespace LinearAlgebra
  {
    LinAlgWorkspaceCpu::LinAlgWorkspaceCpu()
    {
    }

    LinAlgWorkspaceCpu::~LinAlgWorkspaceCpu()
    {
    }

    void LinAlgWorkspaceCpu::initializeHandles()
    {
    }

    void LinAlgWorkspaceCpu::resetLinAlgWorkspace()
    {
      // No resources to reset in CPU workspace
      return;
    }
  } // namespace LinearAlgebra
} // namespace GridKit
