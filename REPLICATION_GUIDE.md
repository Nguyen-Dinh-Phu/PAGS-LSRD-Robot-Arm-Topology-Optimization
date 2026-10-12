# Replication guide

This guide documents the MATLAB source on branch `review/matlab-source-and-replication`. It is not a claim of independently validated numerical replication.

## Prerequisites and branch

Use a MATLAB release supporting `tiledlayout`, `yyaxis`, and `MarkerIndices`. Download the **review branch** (the default branch `Nguyen-Dinh-Phu` may not contain the MATLAB functions):

```bash
git clone --branch review/matlab-source-and-replication --single-branch https://github.com/Nguyen-Dinh-Phu/PAGS-LSRD-Robot-Arm-Topology-Optimization.git
```

Run from MATLAB:

```matlab
cd('path/to/PAGS-LSRD-Robot-Arm-Topology-Optimization')
MAIN
```

`MAIN.m` adds `FUNCTION_PAGS/` to the MATLAB search path.

## Numerical settings

- **Domain:** 120 × 30; **mesh:** 360 × 90 Q4 elements (`Domain(120,30,360,90)`).
- **Loads and supports:** robotic-arm boundary case `'3'`; `ApplyBC(...,'3',5000)` specifies `Fx=0` and total `Fy=-5000 N`. The current routine fixes the left joint's selected inner-rim nodes.
- **Material:** `E11=73.2e9 Pa`, `E22=4.1e9 Pa`, `G12=2.2e9 Pa`, `nu12=0.33`. `Emin=1e-4` is a dimensionless residual stiffness factor.
- **Optimization:** `Vmaxin=0.60`, `Inter=400`, `nvol=200`, `dt=0.1`, `tau=1`, `betaVol=25`, `dVol=-0.001`.
- **History:** baseline `w=0, theta=0°`; four `w=0.4` cases at `theta=[0,30,45,60]°`, with `alpha=1` and `delta=0.6`.
- **Sweep:** 24 combinations of `theta=[0,30,45,60]°` and `w=[0,0.1,0.2,0.3,0.4,0.5]`; `alpha=0.6`, `delta=1`.

Use `study.runMode='history'`, `'sweep'`, or `'all'` to control execution; the uploaded driver defaults to `'all'`.

## Results and display behavior

The code is configured to open the Fig. 3–7 computational plots and the Table 3 topology matrix directly in MATLAB, a total of 10 figure windows when both studies run. No PNG, CSV, or MAT files are written by `PAGS_Output.m`. However, `PAGS_Study.m` **does create an empty `Results/` directory** if missing and prints a message naming that directory. The `studyResults` workspace variable contains the computed result structures.

History arrays record responses evaluated **before** each level-set update, while `Cfinal` and `d2final` refer to the state **after** the final update. The convergence study passes the preceding case's final penalization exponent into the next case's reference-state calculation. The sweep uses `referencePenal=1.0` for each case.

`Calculate_RC0` estimates the reference center; its fallback for non-clamped supports uses a magnitude-weighted reaction centroid rather than the moment-equivalent resultant location. Its computed reference point must be interpreted accordingly.

## Verification before publication

1. Test the orthotropic Q4 element on known displacement/patch cases and check the local stiffness-node ordering against `SYS.edofMat`.
2. Verify consistent element numbering across mesh connectivity, element centroids, passive sets, and density fields.
3. Check force resultant, support DOFs, and selected passive regions numerically.
4. Reconcile the manuscript's published `G12`, `Emin`, material units, loads, and element count with the actual source.
5. Run the full history and 24-case sweep at 360 × 90, compare compliance, CG–RC squared distance, final topology, and reported figures/table to the original research outputs.
6. Document the MATLAB release, compute environment, tested commit SHA, and numerical tolerances.

These steps have **not yet been completed or verified** through execution of the uploaded code.
