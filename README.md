# Adaptation reduces inter-observer metamerism with narrow-band primaries

MATLAB code and data for the simulation study **"Adaptation reduces
inter-observer metamerism with narrow-band primaries."** The pipeline simulates
a large population of individual colorimetric observers, reproduces test
surfaces on a family of virtual RGB displays whose primary bandwidth is swept
from broadband to narrowband, and quantifies how much of the resulting
**observer metamerism** is removed by chromatic adaptation — across three
correlated-colour-temperature (CCT)-matched illuminant classes.

## What the study computes

1. **Individual observers.** 10,000 individual colorimetric observers are drawn
   from the Asano *Individual Colorimetric Observer Model* built on the CIE 2006
   cone-fundamental machinery. Individual differences in lens/ocular-media
   density, macular-pigment density, photopigment optical densities, age, and
   field size follow Asano's continuous distributions. The cone peak
   wavelengths (λmax) are instead **gene-determined and discrete**: the L cone
   is bimodal (Ser180Ala polymorphism, +1.2 nm / −1.5 nm at 0.56 : 0.44), the
   M cone is fixed (monomorphic at codon 180), and the S cone is a symmetric
   two-point distribution (±1.35 nm). The reference is the CIE 1931 2° standard
   colorimetric observer.

2. **Per-observer metameric matches.** For every surface reflectance under every
   test illuminant, each observer's LMS cone excitations to the "world" spectrum
   (illuminant × reflectance) are computed, and the RGB drives that make a
   virtual Gaussian-primary display reproduce *that observer's* LMS are solved.
   A display driven to match one observer is generally not a match for another,
   so the reproduced stimuli disagree once read out in a common space — that
   disagreement is **observer metamerism**.

3. **Metric and adaptation.** Observer metamerism is scored as the CAM16-UCS
   colour difference ΔE′ against the standard observer. Chromatic adaptation is
   modelled with the CIECAM16 / CAT16 von Kries gain at three degrees of
   adaptation D ∈ {0, 0.5, 1}; the observer stimulus and the reference are
   placed at the same D, so the global adaptation-state cast cancels and only
   the between-observer disagreement remains.

The primary bandwidth (σ = 5, 10, 20 nm) and the illuminant class (daylight,
smooth-artificial, peaky-artificial; 12 CCTs each over 4000–7000 K) are the two
controlled variables.

## How to run

Requirements: **MATLAB** (developed with R2025b) and the **Statistics and
Machine Learning Toolbox** (`signrank`, `ranksum`). The **Parallel Computing
Toolbox** is optional and speeds up the main simulation (the code falls back to
serial automatically).

```matlab
main          % full rebuild (~30–35 min): observers → simulation → figures
main(false)   % reuse results/*.mat and only redraw the figures
```

`main` runs from the repository root and drives, in order:
`build_illuminant_sets` → `run_metamerism_analysis` → `compare_nonlinear_adapt_ucs`
→ the six figure scripts. The observer ensemble is generated with a fixed RNG
seed (`rng(0)`), so results are identical from run to run.

## Repository layout

```
main.m                     one-click driver (build results, then draw figures)
src/analysis/              observer generation, the simulation, and statistics
  generate_coneFund_family_AsanoCIE2006.m   individual-observer generator (discrete λmax)
  build_illuminant_sets.m                   daylight / smooth / peaky illuminant classes
  run_metamerism_analysis.m                 the main simulation (writes results/*.mat)
  compare_nonlinear_adapt_ucs.m             pools CAM16-UCS metamerism over the D sweep
  stats_results_narrative.m                 medians and adaptation-reduction statistics
  stats_ellipse_and_daylight.m              between-observer ellipse-area statistics
  fig5_blue_hue_test.m                      blue-surface / short-wavelength test
  optimal_white_upperlimit.m                best-white upper bound on the reduction
  adapting_white_variants.m                 alternative adapting-white anchors (mean LMS, mean chromaticity)
  check_uv_saturation.m                     stimulus-saturation sanity check
src/plotting/              one script per manuscript figure (self-contained)
src/functions/             shared helpers (figure styling, daylight SPD, IES TM-30, …)
data/                      bundled inputs (CIE tables, reflectances, illuminant sets)
results/                   generated .mat outputs   (git-ignored; regenerate with main)
figs/                      rendered figures         (git-ignored; regenerate with main)
```

## Figures (script name ↔ manuscript figure)

The figure **script** numbers were fixed early and no longer match the final
**manuscript** figure order (the figures were renumbered during revision):

| Script                                | Manuscript figure |
|---------------------------------------|-------------------|
| `plot_fig1_observers`                 | Figure 1 — individual-observer cone fundamentals |
| `plot_fig2_stimuli`                   | Figure 2 — reflectances, illuminants, display primaries, gamut |
| `plot_fig3_schematic`                 | Figure 3 — metamerism-quantification schematic |
| `plot_fig4_degree_of_metamerism`      | Figure 4 — metamerism vs degree of adaptation |
| `plot_fig5_ellipses`                  | Figure 5 — CAM16-UCS (a′, b′) ellipses before/after adaptation |
| `plot_fig6_diagonality`               | Figure 6 — cross-observer map Γ distribution |
| `plot_fig7_adapting_white`            | Figure 7 — adapting-white definitions (bar graph, mean LMS / mean chromaticity vs brightest surface) |

Each figure script writes into its own `figs/figN/` subfolder (named by the
*script* number).

## Third-party components

`src/functions/brewermap.m` (colour maps, S. Cobeldick) and
`src/functions/plotxyChromaticity.m` (The MathWorks) retain their original
licence/copyright headers.

## Authors and licence

Takuma Morimoto, Allie C. Schneider, and Hannah E. Smithson.
Released under the MIT Licence (see `LICENSE`).
