# Manuscript outline: Environmental Modelling & Software

Decision (Diogo, 19 Oct 2026): target Environmental Modelling & Software (EMS, Elsevier).
Dissemination to the eDNA community is then ours to do (preprint, talks, direct contact with
the groups running stream additions).

## What EMS expects (from memory of the guide; ScienceDirect blocks automated reading, so
## confirm on the journal's Guide for Authors before submission)

- Article types: research paper (the one for us), short communication, review, position paper.
  Research papers are normally kept under about 8,000 words of main text.
- A mandatory **Software availability** section (sometimes "Software and data availability")
  immediately after the Introduction or before the References, listing: name of software,
  developers, contact address and e-mail, year first available, hardware required, software
  required, programming language, program size, availability (URL, DOI) and cost, licence.
- Highlights: 3-5 bullets, each under 85 characters.
- Abstract about 250 words, no references.
- EMS reviewers weigh: the model is stated completely, the software is open and reproducible,
  the evaluation is honest about limits, and the paper says what the tool is for.
- Elsevier numbered reference style for this journal; hybrid open access (APC if chosen; the
  subscription route costs nothing).

## Title (working)

"A web application for interpreting environmental DNA concentration profiles in streams:
component decay, reactive mobile-immobile transport and source localisation, with a field
evaluation"

Shorter alternative: "An open Shiny application for eDNA removal rates, detection distance and
source location in flowing water"

## Framing

The tool paper the eDNA community does not yet have: four analyses that experimentalists run
after a steady-state addition (removal rates per component, detection distance, hyporheic
versus water-column removal, source distance from composition), in one place, with the models
stated and the code open in two languages. The evaluation is the honest part: the thesis
mesocosm validation (three components, known answer, works) alongside 36 field additions
(two components, known answer by construction, works for rates, partly for source location).
The design message that falls out: measure at least three independently decaying components,
co-release a conservative tracer, and sample densely; with two components the source-location
problem is one-dimensional and its tolerance parameter has no effect.

## Sections

1. **Introduction** (about 900 words). eDNA in flowing water; what existing tools do (eDITH
   for network-scale production, eDNAjoint for occupancy/abundance with traditional surveys,
   the Diana et al. occupancy Shiny app) and what none of them do: use the composition of the
   eDNA signal and the mobile-immobile structure of stream transport. Contribution list.

2. **Models** (about 1,500 words). Exactly the thesis, with the code's conventions stated.
   2.1 Independent first-order decay per component (Eq. 4.5) and why (model comparison in
   Brandao-Dias 2023a, Sec. 4.2). 2.2 Composition and its trajectory (Eqs. 4.6-4.8);
   Algorithm 2 with the implementation rules the code applies: replicate-max normalisation
   before regression, uniform draws of intercept and slope within the 95% CI, slope clipped at
   zero, distance over n-1 fractions, boundary solutions discarded. 2.3 Detection distance.
   2.4 Reactive SMIM at steady state (Eqs. 4.10-4.15), t1, t2 conventions, inversion for r_h.

3. **Software** (about 1,200 words). Architecture of the Shiny app (tabs, inputs, gating by
   data type); input templates; defaults; the Python port (`edna_py/original.py`) as a second
   implementation with the cross-implementation checks (regression, incomplete gamma, k, r_h
   to 12 digits against R); what cannot be matched bit-for-bit (random stream, SANN) and how
   that is handled. Installation, deployment, licence. Table 1: analyses, inputs, outputs.

4. **Evaluation data** (about 700 words). 4.1 Mesocosm validation with three components
   (thesis Sec. 4.3.1.2 and Table 4.2, cited and summarised; the positive control).
   4.2 LEEF 2021: 4 substrates x 4 biofilm phases x 2 species, 10-50 m, two size fractions,
   with rhodamine-derived U, D, Lambda, beta. 4.3 Fort Custer: natural stream, 4 seasons,
   20-100 m, conductivity-corrected. Table 2: datasets. What is missing (bucket compositions)
   and the reference-station design used instead, which is the thesis's own validation design
   transposed to space.

5. **Results** (about 1,800 words).
   5.1 Removal rates: app regression against Snyder's Bayesian hierarchical fits (Fort Custer
   Vdep 0.51/0.80/0.60/0.25 vs 0.53/0.73/0.72/0.33 m/d); size ordering; LEEF significance
   counts; effect of the normalisation convention (median 0.0007 m^-1).
   5.2 Detection distance: LEEF medians 205 m (carp) and 248 m (steelhead) at the ddPCR floor;
   why Fort Custer's conductivity-corrected units make the LOD tab inapplicable there.
   5.3 SMIM: r_s from the Dryad mesocosm series; r_h 2e-3 to 3e-2 1/s for 40 significant
   profiles; r_h/r_s median 190, matching the ARG result; r_h insensitive to r_s, 27%
   sensitive to t2; 21 profiles with k <= 0 cannot be represented; interpretation as
   "removal while in storage".
   5.4 Source localisation (reference-station test): coverage 24/30 LEEF, 3/4 Fort Custer;
   [P10, P90] spans 55-68% of the window; the tolerance has no effect at n = 2 (with proof in
   one paragraph); histogram mode front-loaded; thesis case D comparison (0.58 relative width).

6. **Discussion** (about 1,200 words). What the tool is for and what it is not. Design
   guidance: components, tracer, stations. Why the SMIM r_h absorbs deposition and what an
   explicit deposition term would need. Where the composition method gains information
   (three or more components; mitochondrial gene ratios as an internal reference, thesis p. 98).
   Limitations: single steady source, independent decay, uniform independent parameter draws.
   Roadmap: censored non-detects, joint parameter uncertainty, network extension via eDITH.

7. **Software availability** (mandatory section, half a page).

8. **Conclusions** (about 250 words).

## Figures and tables

- Fig. 1 Concept: component decay, composition trajectory in the simplex, source tolerance
  region (redrawn from thesis Fig. 4.3-4.4; our own drawing).
- Fig. 2 Workflow and screenshots of the four tabs with the Fort Custer winter data loaded.
- Fig. 3 Fort Custer profiles by size with app fits; inset app Vdep vs Snyder Vdep (4 seasons x
  3 components).
- Fig. 4 LEEF r_h by species, size, substrate and phase, with CI bars from the slope CI; panel
  with sensitivity to t2 and r_s.
- Fig. 5 Source-localisation histograms for the four Fort Custer seasons and the coverage
  summary across all 34 datasets; a panel showing why tolerance is irrelevant at n = 2
  (one-dimensional trajectory crossing f_source).
- Fig. 6 Thesis three-component validation summary (Table 4.2 as a figure) next to the
  two-component field result, to make the design message visual.
- Table 1 Analyses, required inputs, outputs, defaults.
- Table 2 Datasets.
- Table 3 r_s by species and fraction from the mesocosm series.
- Table 4 Summary metrics (coverage, relative width, accepted trials) by experiment.
- Supplement: full per-dataset tables (the CSVs in `results/`), the port's cross-checks, the
  app-format input files.

## Open decisions for Diogo

1. Authorship and order. Daniel (app, thesis chapter), Elise Snyder and Jennifer Tank (both
   field datasets), Pedro Brandao-Dias (decay model, mesocosm data), Diogo. Data owners should
   see the field results before anything is written up.
2. Whether the Python port is part of the released software (recommended: yes, as a second
   implementation with the cross-checks in the paper) or stays as analysis code.
3. A versioned release: tag Daniel's GitHub repo, archive on Zenodo with a DOI, confirm the
   shinyapps.io deployment is live and which version it runs. EMS will ask for all three.
4. Licence for the code (the repo has none at present).
5. Whether to run the thesis mesocosm validation again through the port (needs that dataset)
   or cite Table 4.2 as published.

## Writing order

Software availability and Table 1 first (they fix what "the software" is), then Models,
then Results from `07` and `08`, then Software, Data, Discussion, Introduction, Abstract,
Highlights last.
