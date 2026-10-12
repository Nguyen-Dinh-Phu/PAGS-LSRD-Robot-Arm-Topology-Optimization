# Pre-publication review checklist

**Source branch:** `review/matlab-source-and-replication`

**Default branch:** `Nguyen-Dinh-Phu` (MATLAB source not yet merged).

- [x] Upload `MAIN.m` and all 13 MATLAB functions in `FUNCTION_PAGS/` to the review branch.
- [x] Verify that the 14 MATLAB paths exist on GitHub.
- [x] Update README and replication instructions for the actual branch, source layout, and display-only output.
- [ ] Confirm the MATLAB files execute on a specified MATLAB release; check case-sensitive filenames and `addpath`.
- [ ] Validate local Q4 stiffness node order versus `SYS.edofMat` and test mesh element numbering.
- [ ] Confirm the total vertical load `Fy=-5000 N`, fixed DOFs, and non-design element sets.
- [ ] Reconcile `G12` and `Emin` settings with the submitted manuscript.
- [ ] Execute `history` and all 24 `sweep` cases using the 360 × 90 mesh.
- [ ] Compare Fig. 3–7 curves, Table 3 geometry and numerical quantities with original manuscript outputs.
- [ ] Document the verified MATLAB version and commit SHA.
- [ ] Obtain author agreement on release conditions and a suitable software license.
- [ ] Add accurate manuscript citation metadata once final bibliographic details exist.
- [ ] Review changes and merge the verified source into the default branch when ready.

**Current output behavior:** Plot windows are displayed, and `PAGS_Study.m` creates an empty `Results/` folder; `PAGS_Output.m` writes no PNG/CSV/MAT files.

The review branch remains publicly accessible. Completed upload and documentation steps do not constitute numerical validation.
