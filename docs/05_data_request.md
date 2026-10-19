# Data request to Elise Snyder (draft, for Diogo to send)

Everything below is referenced in SNYDERE92024D.pdf but is not in the two data files.

1. Release-solution ("bucket") eDNA concentrations by size fraction (1.2 um and 0.4 um), for each
   addition: LEEF 2021 (4 streams x 4 phases x CARP and STLHD; Sec. 5.3.2 says 50-100 mL bucket
   samples were taken before and after each addition in each stream; Fig. 5.2E-F plots them by
   phase) and Fort Custer (4 seasons; Fig. A.5). Per replicate if available, otherwise the
   before/after pair. These are the "source composition" input for the source-location method.

2. For Fort Custer: plateau and background conductivity at each station and season, so `plteDNA`
   can be converted back to GC/mL and the dilution correction applied explicitly.

3. For LEEF 2021: the truncated-power-law cutoffs t1 and t2 from the rhodamine breakthrough fits
   (Volponi et al., in review), per stream and phase if they vary, alongside the U, D, Lambda, beta
   already in the covariables file.

4. Mesocosm water-column removal rates by size fraction (1.2 and 0.4 um, or 10/1/0.2 um) in units
   of 1/s or 1/h for CARP and STLHD (Snyder et al. 2023; Brandão-Dias et al. 2023a; Ch. 2 and 3),
   to serve as r_s in the SMIM inversion.

5. If easy: the recirculating-mesocosm dataset used in Hallack Ch. 4 Table 4.2 (24 mesocosms, STHD,
   40 min / 6 / 18 / 48 h, 10/1/0.2 um, short and long amplicons; Snyder et al. 2023 /
   Brandão-Dias et al. 2025), so Table 4.2 can be reproduced as a check of the tool.
