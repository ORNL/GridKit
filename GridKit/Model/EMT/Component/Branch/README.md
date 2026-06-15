# Branch Model

## Introduction

EMT branch models represent three-phase network connections between buses in instantaneous abc coordinates.

## Types

### Lumped Parameter

Lumped transmission line models approximate the branch with finite network elements (sometimes referred to as the $\pi$-model). GridKit currently only implements constant parameter.

- `BranchLumpedConstant` (See [BranchLumpedConstant](BranchLumpedConstant/README.md))
- `BranchLumpedFrequencyDependent`

### Distributed Parameter

Distributed transmission line models preserve propagation and delay.

- `BranchDistributedConstant`
- `BranchDistributedFrequencyDependent` (See [BranchDistributedFrequencyDependent](BranchDistributedFrequencyDependent/README.md))
