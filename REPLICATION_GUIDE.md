# Replication guide

## Prerequisites

MATLAB version supporting `tiledlayout`, `MarkerIndices`, and the plotting functions used in `PAGS_Output.m`. No MATLAB execution or full numerical benchmark has yet been performed on the refactored source.

## Execution

```matlab
cd('path/to/PAGS-LSRD-Robot-Arm-Topology-Optimization')
MAIN
```

The entry point adds the project root and `FUNCTION_PAGS` to the MATLAB path. Select the study with `study.runMode`:

- `history`: reference (`w=0`) and `w=0.4` at 0°, 30°, 45°, 60°; `alpha=1.0`, `delta=0.6`.
- `sweep`: 24 combinations of `w=0:0.1:0.5` and `theta=[0,30,45,60]`; `alpha=0.6`, `delta=1.0`.
- `all`: both studies.

The intended iterations per case are `opt.Inter=400`; target volume fraction `opt.Vmaxin=0.60`; continuation length `opt.nvol=200`.

## Fidelity and plotting

**Mesh:** `Domain(120,30,120,30)` in the short driver is a smoke-test mesh, not the published `360×90` discretization. For manuscript-resolution results, use `Domain(120,30,360,90)` and retain the published boundary case and physical parameters.

**Load:** `ApplyBC(...,'3',5000)` is a total resultant `Fy=-5000 N` in the current case-3 code. Earlier input conventions used a different nominal force multiplier: do not mix these.

**Display:** `PAGS_Output` displays figures; it does not save PNG, CSV, or MAT files. To archive scientific numerical outputs, separately store approved result data and explain any changes in the manuscript/repository.

**Sampling:** Iteration histories represent states evaluated **before** each update, while final metrics are evaluated **after** the last update. Expect a possible small difference at iteration 400.

**Interpretation:** `Calculate_RC0` supplies the reference center used for CG–RC comparisons; its branch for non-clamped supports uses a reaction-magnitude-weighted location, not a force/moment-equivalent line of action.

## Verification required before claiming reproduction

1. Check Q4 displacement DOF connectivity versus the local `coord` ordering used to integrate `SYS.KE`. The original assembly may have inconsistent node ordering and should not be silently changed.
2. Confirm `mesh.info.Centroid`, the passive sets and `SYS.edofMat` share one element enumeration.
3. Verify resultant force, fixed DOFs, and passive element maps for boundary case 3.
4. Compare baseline compliance, volume history, all four angles, and 24-case sweep results against the original MATLAB scripts and manuscript tables at the **same mesh**.
5. Record MATLAB version, compute environment and tested commit SHA.

Until these checks have passed, the repository should be described as a **research implementation under numerical verification**, not as independently validated exact replication.
