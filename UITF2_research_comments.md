# UITF2 Research Comments and Validation Notes

## Study Rationale

This v2 UITF simulation package supports a DNP target-material irradiation
study. The central physics question is whether an 8 MeV electron beam at UITF
can deliver a useful, spatially interpretable dose to cryogenic ND3/NH3 target
material while preserving a practical heat-load margin. The target must be
irradiated from both sides because the electron range and depth-dose falloff
make a one-sided exposure an incomplete representation of the proposal
protocol.

The simulation should therefore be judged on two levels:

1. Qualitative validity: the beam, geometry, front/back logic, spectra, and
   dose shapes must look physically plausible.
2. Quantitative validity: normalization, histories, fluence, dose, heat load,
   and target mass must be traceable from TOPAS scorer output to reported
   values.

## Production Definition

- Front side: `UITF2-10M-front.txt`
- Back side: `UITF2-10M-back.txt`
- Histories: 10,000,000 per side
- Front output: `UITF_DNP_Output_10M_front`
- Back output: `UITF_DNP_Output_10M_back`
- Combined interpretation: front run plus target-coordinate-oriented back run
- Physical timing: front and back are sequential irradiations, not simultaneous

Legacy files named `UITF2-100M-front.txt` and `UITF2-100M-back.txt` remain in
the directory for compatibility, but they now run the v2 10M-per-side protocol.

## Research Comments for Manuscript Reporting

The simulation should be written up as a validation chain, not as a figure
gallery. A strong manuscript narrative is:

1. Establish the experimental/proposal geometry: beam energy, raster area,
   cryostat, target length, target mass, and material state.
2. Show why a one-sided run is insufficient by presenting the front-only
   depth-dose falloff.
3. Add the back-side run in target coordinates and show the cumulative
   front+back dose.
4. Validate beam transport using entrance/exit fluence and entrance spectrum.
5. Validate spatial coverage using the R-Z dose map and radial dose scorer.
6. Convert per-primary scorer values to physical dose and heat load.
7. State which quantities are directly scored and which are derived proxies.
8. Close with limitations and the measurements needed for benchmarking.

## Qualitative Validity Checklist

- The front dose profile has a believable depth falloff for 8 MeV electrons.
- The back profile complements the front profile after orientation.
- The combined profile is smoother and more proposal-relevant than either
  one-sided profile.
- The 2D dose map shows a rastered volume, not a narrow un-rastered beam.
- The entrance energy spectrum remains consistent with the nominal 8 MeV beam.
- Bremsstrahlung/gamma outputs are treated as secondary radiation context.
- The report never implies simultaneous front/back heat load.

## Quantitative Validity Checklist

- Each production run uses exactly 10,000,000 histories.
- TOPAS scorer sums are divided by 10,000,000 before physical scaling.
- Beam current is converted with the exact electron charge.
- Dose-rate and cumulative-dose calculations use the correct per-side
  irradiation times.
- Heat load is derived from total deposited energy whenever possible and is
  cross-checked against dose-to-mass conversion.
- Target mass is reported as 250 mg, with the distinction between proposal
  bulk density and solid-density microdosimetry made explicit.
- Any radical-yield or LET-like conclusion is identified as a derived model
  proxy, not a direct measurement.

## Suggested Figure Order

1. Geometry/protocol schematic or annotated TOPAS geometry screenshot.
2. Front-only and back-only depth-dose profiles.
3. Front+back cumulative depth-dose profile.
4. R-Z dose map for cumulative front+back exposure.
5. Entrance and exit fluence/spectral validation.
6. Heat load versus current with cold and warm operating points.
7. Derived radical-production proxy with clear caveats.
8. Summary table of histories, normalization, fluence, dose, and heat load.

## Minimum Claim Discipline

Use strong language for direct Monte Carlo quantities:

- dose per primary;
- energy deposited per primary;
- electron/photon fluence;
- front/back cumulative dose after superposition.

Use cautious language for inferred material-response quantities:

- radical yield;
- polarization performance;
- low-temperature radical stability;
- LET-weighted chemical effectiveness.

Those inferred quantities become formal claims only after comparison with an
experimental benchmark or a clearly cited empirical model.
