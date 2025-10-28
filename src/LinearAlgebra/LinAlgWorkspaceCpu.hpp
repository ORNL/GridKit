#pragma once

namespace GridKit
{
  namespace LinearAlgebra
  {
    class LinAlgWorkspaceCpu
    {
    public:
      LinAlgWorkspaceCpu();
      ~LinAlgWorkspaceCpu();
      void initializeHandles();
      void resetLinAlgWorkspace();
    };
  } // namespace LinearAlgebra
} // namespace GridKit
