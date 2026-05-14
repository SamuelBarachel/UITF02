# UITF02 DNP Target Irradiation Study

This repository contains a TOPAS/Geant4 Monte Carlo study of two-sided
8 MeV electron irradiation for cryogenic packed ND3 target-material
preparation at the Jefferson Lab Upgraded Injector Test Facility (UITF).

The central question is whether a front/back irradiation protocol can deliver
a traceable target-volume dose field while keeping the cold-operation heat load
in a plausible range. The analysis reports quantities directly supported by
Monte Carlo transport: dose, fluence, spectral diagnostics, energy deposition,
radial coverage, and target heat load.

## Report

- `Project_Report.md` is the main report source.
- `pandoc_uitf2.yaml` and `uitf2_manuscript.css` support HTML rendering.
- `outputs/manual-20260507-uitf/presentations/uitf-dnp-topas/output/UITF_DNP_TOPAS_MC_Report.pptx`
  is the presentation deck prepared from the analysis figures.

The report includes a University of New Hampshire title page, abstract,
methods, results and analysis, discussion, conclusion, acknowledgments, data
availability statement, and references.

## Key Results

- The 6.5 cm target length is not well represented by one-sided 8 MeV electron
  exposure alone.
- Sequential front/back irradiation gives the proposal-relevant cumulative
  target-coordinate dose field.
- The cold 1 uA case delivers approximately `4.97e15 e-/cm2` per side, with
  estimated heat load `0.064 W` per active side.
- The warm 10 uA case delivers approximately `9.97e16 e-/cm2` per side, with
  estimated heat load `0.644 W` per active side.
- Radical chemistry and final DNP polarization are not directly modeled by
  TOPAS and require ESR, NMR, or polarization benchmarking.

## Main Files

- `UITF2-10M-front.txt` and `UITF2-10M-back.txt`: production TOPAS cards for
  the front and back irradiation directions.
- `UITF2-050426.txt`: main TOPAS geometry and scorer definition.
- `UITF2G-050426.txt`: OpenGL geometry viewer wrapper.
- `analyze_DNP_UITF2.py`: single-run scorer post-processing.
- `analyze_DNP_UITF2_MATLAB.m`: combined front/back analysis and report figure
  generation.
- `figures_matlab/pptx_deck_assets/`: final figure assets used by the report
  and presentation.

## Figure Preview

![Front, back, and combined cold cumulative dose profiles](figures_matlab/pptx_deck_assets/01_front_back_superimposed_cold_dose.png)

## Reproducing the Analysis

Run the front and back production cards, then run the post-processing scripts:

```sh
./run_UITF2.command
```

The wrapper expects TOPAS and MATLAB to be installed in the application paths
defined near the top of `run_UITF2.command`. The Python analysis requires
`numpy`, `pandas`, `matplotlib`, and `scipy`.

## Repository Notes

Local collaboration documents and external reference PDFs are intentionally not
tracked in this public repository. The cited literature remains listed in the
report references.
