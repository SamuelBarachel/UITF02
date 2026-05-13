%% analyze_DNP_UITF2_MATLAB
% MATLAB post-processing for the UITF2 TOPAS / Geant4 DNP target simulation.
%
% This version is front/back aware. It loads:
%   UITF_DNP_Output_10M_front
%   UITF_DNP_Output_10M_back
%
% If both folders exist, the analysis is already superimposed:
%   - front-side dose is plotted by itself
%   - back-side dose is plotted by itself
%   - front + back cumulative dose is plotted as the combined proposal result
%
% If those folders do not exist, it falls back to the current folder so debug
% CSV files can still be inspected.

clear; clc; close all;

%% Configuration
cfg = struct();
cfg.eChargeC = 1.602176634e-19;
cfg.mevToJ = 1.602176634e-13;
cfg.beamEnergyMeV = 8.0;
cfg.nPrimaries = 10000000;        % v2 production histories per irradiation side

cfg.targetLengthCm = 6.5;
cfg.targetRadiusCm = 1.42;
cfg.interactionAreaCm2 = 6.42;
cfg.targetMassMg = 250.0;

cfg.warmCurrentUA = 10.0;
cfg.coldCurrentUA = 1.0;
cfg.warmFluenceTarget = 1.0e17;   % e-/cm2 per side
cfg.coldFluenceTarget = 5.0e15;   % e-/cm2 per side
cfg.warmHoursPerSide = 2.85;
cfg.coldHoursPerSide = 1.42;
cfg.solidND3Density = 1.007;      % g/cm3, for CSDA estimate
cfg.bulkDensity = (cfg.targetMassMg * 1e-3) / ...
    (pi * cfg.targetRadiusCm^2 * cfg.targetLengthCm);
cfg.avogadro = 6.02214076e23;
cfg.nd3MolarMassGmol = 14.0067 + 3 * 2.014101778;
cfg.nd3ElectronsPerMolecule = 10;
cfg.packedND3ElectronDensityCm3 = electronDensityND3(cfg.bulkDensity, cfg);
cfg.solidND3ElectronDensityCm3 = electronDensityND3(cfg.solidND3Density, cfg);

scriptDir = string(fileparts(mfilename('fullpath')));
if strlength(scriptDir) == 0, scriptDir = string(pwd); end
frontDir = fullfile(scriptDir, "UITF_DNP_Output_10M_front");
backDir = fullfile(scriptDir, "UITF_DNP_Output_10M_back");
figDir = fullfile(scriptDir, "figures_matlab");
if ~exist(figDir, "dir"), mkdir(figDir); end
singleFigDir = fullfile(figDir, "pptx_deck_assets");
if ~exist(singleFigDir, "dir"), mkdir(singleFigDir); end

% Bright, high-contrast presentation palette. Every text and axes color is set
% explicitly so MATLAB dark mode cannot make titles or labels disappear.
pal = struct();
pal.page = [0.965 0.985 1.000];        % pale blue page
pal.panel = [1.000 0.995 0.955];       % warm light plotting panel
pal.ink = [0.050 0.090 0.220];         % navy text, never black
pal.coldFront = [0.000 0.325 0.760];   % saturated blue
pal.coldBack = [0.000 0.650 0.420];    % saturated green
pal.warmFront = [0.920 0.330 0.000];   % saturated orange
pal.warmBack = [0.820 0.180 0.600];    % magenta
pal.combined = [0.500 0.100 0.920];    % bright purple combined curve
pal.reference = [0.360 0.300 0.620];   % muted violet reference lines
pal.barBlue = [0.000 0.450 0.850];
pal.barOrange = [0.930 0.420 0.050];
pal.barGreen = [0.000 0.620 0.420];
pal.barPink = [0.820 0.180 0.600];
pal.edge = [0.070 0.120 0.280];

haveFront = isfolder(frontDir) && isfile(fullfile(frontDir, "DoseVsDepth_ND3_2K.csv"));
haveBack = isfolder(backDir) && isfile(fullfile(backDir, "DoseVsDepth_ND3_2K.csv"));

if haveFront
    front = loadTopasRun(frontDir, "Front", cfg);
else
    warning("Front output folder not found. Falling back to current folder.");
    front = loadTopasRun(".", "SingleRun", cfg);
end

if haveBack
    back = loadTopasRun(backDir, "Back", cfg);
else
    back = [];
end

zMm = front.zMm;
[csdaRangeCm, csdaRangeGcm2] = csdaRangeEstimate(cfg.beamEnergyMeV, cfg.solidND3Density);

frontCold = absoluteDose(front.doseDepthGyPerPrimary, cfg.coldCurrentUA, "cold", 1, cfg);
frontWarm = absoluteDose(front.doseDepthGyPerPrimary, cfg.warmCurrentUA, "warm", 1, cfg);

if ~isempty(back)
    [backDoseForCombination, backMirrorNote] = orientBackDoseForCombination( ...
        front.doseDepthGyPerPrimary, back.doseDepthGyPerPrimary);
    backCold = absoluteDose(backDoseForCombination, cfg.coldCurrentUA, "cold", 1, cfg);
    backWarm = absoluteDose(backDoseForCombination, cfg.warmCurrentUA, "warm", 1, cfg);
    combinedDoseGyPerPrimary = front.doseDepthGyPerPrimary + backDoseForCombination;
    combinedCold = absoluteDose(combinedDoseGyPerPrimary, cfg.coldCurrentUA, "cold", 1, cfg);
    combinedWarm = absoluteDose(combinedDoseGyPerPrimary, cfg.warmCurrentUA, "warm", 1, cfg);
    isSuperimposed = true;
    isEstimatedBack = false;
else
    backDoseForCombination = flipud(front.doseDepthGyPerPrimary);
    backMirrorNote = "No back run loaded; using mirrored front profile as an estimated back-side dose.";
    backCold = absoluteDose(backDoseForCombination, cfg.coldCurrentUA, "cold", 1, cfg);
    backWarm = absoluteDose(backDoseForCombination, cfg.warmCurrentUA, "warm", 1, cfg);
    combinedDoseGyPerPrimary = front.doseDepthGyPerPrimary + backDoseForCombination;
    combinedCold = absoluteDose(combinedDoseGyPerPrimary, cfg.coldCurrentUA, "cold", 1, cfg);
    combinedWarm = absoluteDose(combinedDoseGyPerPrimary, cfg.warmCurrentUA, "warm", 1, cfg);
    isSuperimposed = true;
    isEstimatedBack = true;
end

[combinedDoseMapGyPerPrimary, doseMapNote] = combineDoseMaps(front, back);
frontElectronFluence = depthVector(front.fluenceDepth, front.doseDepthGyPerPrimary);
if ~isempty(frontElectronFluence)
    if ~isempty(back)
        backElectronFluence = depthVector(back.fluenceDepth, front.doseDepthGyPerPrimary);
        if ~isempty(backElectronFluence)
            [backElectronFluence, ~] = orientBackDoseForCombination( ...
                frontElectronFluence, backElectronFluence);
            combinedElectronFluence = frontElectronFluence + backElectronFluence;
            electronFluenceNote = "Front electron fluence + back-oriented electron fluence";
        else
            combinedElectronFluence = frontElectronFluence;
            electronFluenceNote = "Back electron fluence scorer missing; front run only";
        end
    else
        combinedElectronFluence = frontElectronFluence + flipud(frontElectronFluence);
        electronFluenceNote = "Front electron fluence + mirrored-front estimate";
    end
else
    combinedElectronFluence = [];
    electronFluenceNote = "Electron fluence scorer missing";
end

warmFluencePerSide = deliveredFluence(cfg.warmCurrentUA, "warm", 1, cfg);
coldFluencePerSide = deliveredFluence(cfg.coldCurrentUA, "cold", 1, cfg);
warmFluenceBoth = deliveredFluence(cfg.warmCurrentUA, "warm", 2, cfg);
coldFluenceBoth = deliveredFluence(cfg.coldCurrentUA, "cold", 2, cfg);

frontColdHeatW = heatLoadFromEDep(front.totalEdepMeVPerPrimary, cfg.coldCurrentUA, cfg);
frontWarmHeatW = heatLoadFromEDep(front.totalEdepMeVPerPrimary, cfg.warmCurrentUA, cfg);
if ~isempty(back)
    backColdHeatW = heatLoadFromEDep(back.totalEdepMeVPerPrimary, cfg.coldCurrentUA, cfg);
    backWarmHeatW = heatLoadFromEDep(back.totalEdepMeVPerPrimary, cfg.warmCurrentUA, cfg);
else
    backColdHeatW = frontColdHeatW;
    backWarmHeatW = frontWarmHeatW;
end

%% Console summary
fprintf("\n========================================================================\n");
fprintf("  MATLAB SUMMARY - UITF2 FRONT/BACK TOPAS ANALYSIS\n");
fprintf("========================================================================\n");
fprintf("  Front folder:             %s\n", front.sourceDir);
if ~isempty(back)
    fprintf("  Back folder:              %s\n", back.sourceDir);
else
    fprintf("  Back folder:              not loaded; mirrored-front estimate used\n");
end
fprintf("  Front/back superimposed:  %s\n", yesNo(isSuperimposed));
fprintf("  Back-dose orientation:    %s\n", backMirrorNote);
fprintf("  Beam energy:              %.1f MeV electrons\n", cfg.beamEnergyMeV);
fprintf("  Target geometry:          L = %.2f cm, R = %.2f cm\n", cfg.targetLengthCm, cfg.targetRadiusCm);
fprintf("  Interaction area:         %.2f cm2 (proposal: 6.42 cm2)\n", cfg.interactionAreaCm2);
fprintf("  Proposal batch mass:      %.0f mg\n", cfg.targetMassMg);
fprintf("  Bulk-average density:     %.5f g/cm3\n", cfg.bulkDensity);
fprintf("  Packed ND3 electron density: %.3e e-/cm3\n", cfg.packedND3ElectronDensityCm3);
fprintf("  CSDA range in solid ND3:  %.1f mm (%.2f g/cm2)\n", csdaRangeCm * 10, csdaRangeGcm2);
fprintf("\n  Fluence agreement:\n");
fprintf("    Warm proposal/side:     %.2e e-/cm2\n", cfg.warmFluenceTarget);
fprintf("    Warm delivered/side:    %.2e e-/cm2 (%.2fx target)\n", warmFluencePerSide, warmFluencePerSide / cfg.warmFluenceTarget);
fprintf("    Warm front+back total:  %.2e e-/cm2\n", warmFluenceBoth);
fprintf("    Cold proposal/side:     %.2e e-/cm2\n", cfg.coldFluenceTarget);
fprintf("    Cold delivered/side:    %.2e e-/cm2 (%.2fx target)\n", coldFluencePerSide, coldFluencePerSide / cfg.coldFluenceTarget);
fprintf("    Cold front+back total:  %.2e e-/cm2\n", coldFluenceBoth);
fprintf("\n  Dose interpretation:\n");
fprintf("    Front and back are NOT physically simultaneous.\n");
fprintf("    Combined dose = front profile + mirrored/back-oriented profile.\n");
fprintf("    A valid combined curve should not have the same shape as one individual side.\n");
fprintf("    Heat load is quoted per irradiation side, not front+back at the same time.\n");
fprintf("    Electron fluence and material electron density are plotted without normalization.\n");
fprintf("\n  Cold, 1 uA:\n");
fprintf("    Front mean cumulative dose:    %.3e Gy\n", mean(frontCold.cumulativeGy, "omitnan"));
if ~isempty(back) || isEstimatedBack
    fprintf("    Back mean cumulative dose:     %.3e Gy\n", mean(backCold.cumulativeGy, "omitnan"));
end
fprintf("    Combined mean cumulative dose: %.3e Gy\n", mean(combinedCold.cumulativeGy, "omitnan"));
fprintf("    Front heat load:               %.4f W\n", frontColdHeatW);
if ~isempty(back), fprintf("    Back heat load:                %.4f W\n", backColdHeatW); end
fprintf("\n  Warm, 10 uA:\n");
fprintf("    Front mean cumulative dose:    %.3e Gy\n", mean(frontWarm.cumulativeGy, "omitnan"));
if ~isempty(back) || isEstimatedBack
    fprintf("    Back mean cumulative dose:     %.3e Gy\n", mean(backWarm.cumulativeGy, "omitnan"));
end
fprintf("    Combined mean cumulative dose: %.3e Gy\n", mean(combinedWarm.cumulativeGy, "omitnan"));
fprintf("    Front heat load:               %.4f W\n", frontWarmHeatW);
if ~isempty(back), fprintf("    Back heat load:                %.4f W\n", backWarmHeatW); end
fprintf("========================================================================\n\n");

%% Export summary table
backColdMean = mean(backCold.cumulativeGy, "omitnan");
backWarmMean = mean(backWarm.cumulativeGy, "omitnan");

summaryTable = table( ...
    ["Cold"; "Warm"], ...
    [cfg.coldCurrentUA; cfg.warmCurrentUA], ...
    [cfg.coldFluenceTarget; cfg.warmFluenceTarget], ...
    [coldFluencePerSide; warmFluencePerSide], ...
    [coldFluenceBoth; warmFluenceBoth], ...
    [mean(frontCold.cumulativeGy, "omitnan"); mean(frontWarm.cumulativeGy, "omitnan")], ...
    [backColdMean; backWarmMean], ...
    [mean(combinedCold.cumulativeGy, "omitnan"); mean(combinedWarm.cumulativeGy, "omitnan")], ...
    [frontColdHeatW; frontWarmHeatW], ...
    [backColdHeatW; backWarmHeatW], ...
    'VariableNames', ["Mode", "Current_uA", "ProposalFluencePerSide_e_cm2", ...
    "DeliveredFluencePerSide_e_cm2", "DeliveredFluenceFrontBackTotal_e_cm2", ...
    "FrontMeanCumulativeDose_Gy", "BackMeanCumulativeDose_Gy", ...
    "CombinedMeanCumulativeDose_Gy", "FrontHeatLoad_W", "BackHeatLoad_W"]);
summaryTextPath = fullfile(figDir, "UITF2_front_back_MATLAB_summary.txt");
writetable(summaryTable, summaryTextPath, "Delimiter", "\t");

%% Legacy overview export disabled
% The report now uses the single tabbed browser below. The old multi-panel
% overview figures were removed from the interactive workflow so MATLAB does
% not open a stack of separate windows.
fprintf("MATLAB summary table written to: %s\n", summaryTextPath);

%% PPTX-ready tabbed report viewer
% One MATLAB window is created with clickable tabs. Each tab is also exported
% as a separate PNG for the later PPTX report. Front/back superposition is
% used for dose-depth, electron fluence, and the 2D dose map.

figTabs = figure("Color", pal.page, "InvertHardcopy", "off", ...
    "Name", "UITF2 PPTX Figure Browser", "NumberTitle", "off", ...
    "Position", [80 80 1180 760]);
tabGroup = uitabgroup(figTabs);

% 1. Explicit front/back superposition tab.
axS = nextPlotTab(tabGroup, "01 Superimposed", pal);
plot(axS, zMm, frontCold.cumulativeGy, "LineWidth", 3.2, "Color", pal.coldFront); hold(axS, "on");
plot(axS, zMm, backCold.cumulativeGy, "--", "LineWidth", 3.2, "Color", pal.coldBack);
plot(axS, zMm, combinedCold.cumulativeGy, "-.", "LineWidth", 4.0, "Color", pal.combined);
xline(axS, csdaRangeCm * 10, "--", "CSDA", "Color", pal.reference, "LineWidth", 1.8);
xline(axS, cfg.targetLengthCm * 5, ":", "Flip midpoint", "Color", pal.reference, "LineWidth", 1.8);
xlabel(axS, "Position in target from front face (mm)", "Color", pal.ink);
ylabel(axS, "Cold cumulative dose (Gy)", "Color", pal.ink);
title(axS, "Front + Back Dose Superposition", "Color", pal.combined, "FontWeight", "bold");
styleAxes(axS, pal);
if isempty(back)
    lg = legend(axS, "Front run", "Mirrored front estimate", ...
        "Front + mirrored estimate", "Location", "best");
else
    lg = legend(axS, "Front run", "Back run", "Front + back combined", "Location", "best");
end
styleLegend(lg, pal);
saveAxesCopy(axS, singleFigDir, "01_front_back_superimposed_cold_dose");

% 2. Combined cold dose-depth profile tab.
axS = nextPlotTab(tabGroup, "02 Cold Dose", pal);
plot(axS, zMm, combinedCold.cumulativeGy, "LineWidth", 4.0, "Color", pal.combined); hold(axS, "on");
xline(axS, csdaRangeCm * 10, "--", "CSDA", "Color", pal.reference, "LineWidth", 1.8);
xline(axS, cfg.targetLengthCm * 5, ":", "Flip midpoint", "Color", pal.reference, "LineWidth", 1.8);
xlabel(axS, "Depth in target (mm)", "Color", pal.ink);
ylabel(axS, "Cold cumulative dose, front + back (Gy)", "Color", pal.ink);
title(axS, "Combined Cold Dose-Depth Profile", "Color", pal.combined, "FontWeight", "bold");
styleAxes(axS, pal);
saveAxesCopy(axS, singleFigDir, "02_combined_cold_dose_depth");

% 3. Superimposed 2D dose map tab, if available.
if ~isempty(combinedDoseMapGyPerPrimary)
    axS = nextPlotTab(tabGroup, "03 Dose Map", pal);
    imagesc(axS, zMm, front.rMm, combinedDoseMapGyPerPrimary);
    set(axS, "YDir", "normal");
    colormap(axS, turbo);
    cb = colorbar(axS);
    cb.Label.String = "DoseToMedium, front + back (Gy / primary)";
    cb.Label.Color = pal.ink;
    cb.Color = pal.ink;
    xlabel(axS, "Depth z (mm)", "Color", pal.ink);
    ylabel(axS, "Radius r (mm)", "Color", pal.ink);
    title(axS, "2D Dose Map: Front + Back Superposition", "Color", pal.warmFront, "FontWeight", "bold");
    styleAxes(axS, pal);
    text(axS, 0.02, 0.96, doseMapNote, "Units", "normalized", ...
        "Color", pal.ink, "FontWeight", "bold", "BackgroundColor", pal.panel);
    saveAxesCopy(axS, singleFigDir, "03_2d_dose_map_superimposed");
end

% 4. Electron fluence depth tab.
if ~isempty(combinedElectronFluence)
    axS = nextPlotTab(tabGroup, "04 e- Fluence", pal);
    plot(axS, zMm, combinedElectronFluence, "LineWidth", 3.4, "Color", pal.barBlue);
    xlabel(axS, "Depth in target (mm)", "Color", pal.ink);
    ylabel(axS, "Electron fluence (/mm2 per primary)", "Color", pal.ink);
    title(axS, "Electron Fluence Depth Profile, Not Normalized", "Color", pal.barBlue, "FontWeight", "bold");
    styleAxes(axS, pal);
    text(axS, 0.02, 0.96, electronFluenceNote, "Units", "normalized", ...
        "Color", pal.ink, "FontWeight", "bold", "BackgroundColor", pal.panel);
    saveAxesCopy(axS, singleFigDir, "04_electron_fluence_depth_unnormalized");
end

% 5. Material electron density tab.
electronDensityPacked = repmat(cfg.packedND3ElectronDensityCm3, size(zMm));
electronDensitySolid = repmat(cfg.solidND3ElectronDensityCm3, size(zMm));
axS = nextPlotTab(tabGroup, "05 e- Density", pal);
semilogy(axS, zMm, electronDensityPacked, "LineWidth", 3.6, "Color", pal.combined); hold(axS, "on");
semilogy(axS, zMm, electronDensitySolid, "--", "LineWidth", 2.8, "Color", pal.coldBack);
xlabel(axS, "Depth in target (mm)", "Color", pal.ink);
ylabel(axS, "Material electron density (e-/cm3)", "Color", pal.ink);
title(axS, "ND3 Electron Density Used by the Simulation", "Color", pal.combined, "FontWeight", "bold");
styleAxes(axS, pal);
lg = legend(axS, "Packed ND3 target volume", "Solid ND3 material reference", "Location", "best");
styleLegend(lg, pal);
text(axS, 0.02, 0.18, "Computed from ND3 composition, density, and Avogadro's number; not normalized.", ...
    "Units", "normalized", "Color", pal.ink, "FontWeight", "bold", "BackgroundColor", pal.panel);
saveAxesCopy(axS, singleFigDir, "05_material_electron_density_unnormalized");

% 6. Lateral dose uniformity tab using the raw scorer values.
if ~isempty(front.radialDoseGyPerPrimary)
    axS = nextPlotTab(tabGroup, "06 Uniformity", pal);
    radialMean = mean(front.radialDoseGyPerPrimary, "omitnan");
    plot(axS, front.radialMm, front.radialDoseGyPerPrimary, "o-", "LineWidth", 3.2, ...
        "MarkerSize", 7, "Color", pal.coldFront, "MarkerFaceColor", pal.coldFront); hold(axS, "on");
    yline(axS, radialMean * 1.05, "--", "+5% of mean", "Color", pal.warmFront, "LineWidth", 1.6);
    yline(axS, radialMean * 0.95, "--", "-5% of mean", "Color", pal.warmFront, "LineWidth", 1.6);
    xlabel(axS, "Radial position r (mm)", "Color", pal.ink);
    ylabel(axS, "DoseToMedium (Gy / primary)", "Color", pal.ink);
    title(axS, "Raster Coverage / Lateral Dose Uniformity", "Color", pal.coldFront, "FontWeight", "bold");
    styleAxes(axS, pal);
    lg = legend(axS, "Raw radial dose", "+5% of mean", "-5% of mean", "Location", "best");
    styleLegend(lg, pal);
    saveAxesCopy(axS, singleFigDir, "06_lateral_dose_uniformity_raw");
end

% 7. Electron and photon spectra tab.
axS = nextPlotTab(tabGroup, "07 Spectra", pal);
eCenters = linspace(0, 8.1, 160);
plotSpectrum(axS, eCenters, front.energySpecIn, pal.coldFront); hold(axS, "on");
plotSpectrum(axS, eCenters, front.bremsSpecExit, pal.warmFront);
xline(axS, cfg.beamEnergyMeV, ":", "8 MeV nominal", "Color", pal.reference, "LineWidth", 1.8);
xlabel(axS, "Energy (MeV)", "Color", pal.ink);
ylabel(axS, "Energy fluence scorer sum", "Color", pal.ink);
title(axS, "Electron Entrance and Bremsstrahlung Spectra", "Color", pal.combined, "FontWeight", "bold");
styleAxes(axS, pal);
lg = legend(axS, "Entrance electrons", "Exit photons", "Nominal beam energy", "Location", "best");
styleLegend(lg, pal);
saveAxesCopy(axS, singleFigDir, "07_energy_spectra");

% 8. Heat load versus beam current tab.
axS = nextPlotTab(tabGroup, "08 Heat", pal);
currentsUA = logspace(-1, 1, 80);
heatCurveW = arrayfun(@(I) heatLoadFromEDep(front.totalEdepMeVPerPrimary, I, cfg), currentsUA);
loglog(axS, currentsUA, heatCurveW, "LineWidth", 3.4, "Color", pal.warmFront); hold(axS, "on");
yline(axS, 0.1, "--", "0.1 W", "Color", pal.coldFront, "LineWidth", 1.8);
yline(axS, 0.5, "--", "0.5 W", "Color", pal.warmBack, "LineWidth", 1.8);
xline(axS, cfg.coldCurrentUA, ":", "Cold 1 uA", "Color", pal.coldFront, "LineWidth", 1.8);
xline(axS, cfg.warmCurrentUA, ":", "Warm 10 uA", "Color", pal.warmFront, "LineWidth", 1.8);
xlabel(axS, "Beam current (uA)", "Color", pal.ink);
ylabel(axS, "Target heat load per side (W)", "Color", pal.ink);
title(axS, "Target Heat Load vs. Beam Current", "Color", pal.warmFront, "FontWeight", "bold");
styleAxes(axS, pal);
saveAxesCopy(axS, singleFigDir, "08_heat_load_vs_current");

% 9. Acronym and interpretation notes tab.
axS = nextPlotTab(tabGroup, "09 Acronyms", pal);
axis(axS, "off");
notes = [
    "Acronyms and terms used in this analysis"
    ""
    "TOPAS: Tool for Particle Simulation, the Geant4-based Monte Carlo interface used here."
    "Geant4: particle-transport simulation toolkit underneath TOPAS."
    "MC: Monte Carlo, a particle-by-particle stochastic simulation method."
    "UITF: Upgraded Injector Test Facility."
    "DNP: Dynamic Nuclear Polarization, the target-polarization technique."
    "ND3: deuterated ammonia target material."
    "e-: electron. Electron fluence is the TOPAS scorer per primary; electron density is material e-/cm3."
    "CSDA: Continuous Slowing Down Approximation, used for the approximate electron range marker."
    "Gy: gray, absorbed dose in joule per kilogram."
    "MeV: mega-electron-volt, particle energy unit."
    "uA: microampere, beam current."
    "R-Z: cylindrical radius-depth map."
    "LHe / LAr: liquid helium / liquid argon cooling baths."
    "PPTX: PowerPoint report format."
];
text(axS, 0.03, 0.97, notes, "Units", "normalized", "VerticalAlignment", "top", ...
    "Color", pal.ink, "FontWeight", "bold", "FontSize", 12, "Interpreter", "none");
saveAxesCopy(axS, singleFigDir, "09_acronyms_and_interpretation_notes");

savefig(figTabs, fullfile(singleFigDir, "UITF2_tabbed_pptx_figure_browser.fig"));
fprintf("PPTX tabbed figure browser written to: %s\n", fullfile(singleFigDir, "UITF2_tabbed_pptx_figure_browser.fig"));
fprintf("PPTX PNG slide assets written to: %s\n", singleFigDir);

%% Local helper functions
function run = loadTopasRun(sourceDir, label, cfg)
    doseDepthSum = scorerValues(fullfile(sourceDir, "DoseVsDepth_ND3_2K.csv"));
    if isempty(doseDepthSum)
        error("No dose-depth data found in %s", sourceDir);
    end
    run = struct();
    run.sourceDir = sourceDir;
    run.label = label;
    run.doseDepthGyPerPrimary = doseDepthSum(:) ./ cfg.nPrimaries;
    nZ = numel(run.doseDepthGyPerPrimary);
    run.zMm = linspace(0.5, cfg.targetLengthCm * 10 - 0.5, nZ).';
    run.energySpecIn = scorerValues(fullfile(sourceDir, "EnergySpectrum_Entrance.csv"));
    run.bremsSpecExit = scorerValues(fullfile(sourceDir, "BremsSpectrum_Exit.csv"));
    run.edepDepth = scorerValues(fullfile(sourceDir, "EnergyDeposit_DepthProfile.csv")) ./ cfg.nPrimaries;
    run.fluenceDepth = scorerValues(fullfile(sourceDir, "Fluence_DepthProfile.csv")) ./ cfg.nPrimaries;
    run.doseMapGyPerPrimary = doseMapFromCsv(fullfile(sourceDir, "DoseMap_RZ_ND3.csv"), cfg);
    if ~isempty(run.doseMapGyPerPrimary)
        run.rMm = linspace(0.5, cfg.targetRadiusCm * 10 - 0.5, size(run.doseMapGyPerPrimary, 1));
    else
        run.rMm = [];
    end
    radialDose = scorerValues(fullfile(sourceDir, "DoseVsRadius_ND3.csv")) ./ cfg.nPrimaries;
    if ~isempty(radialDose)
        run.radialDoseGyPerPrimary = radialDose(:);
        run.radialDoseNorm = normalizeVector(radialDose(:));
        run.radialMm = linspace(0.5, cfg.targetRadiusCm * 10 - 0.5, numel(radialDose));
    else
        run.radialDoseGyPerPrimary = [];
        run.radialDoseNorm = [];
        run.radialMm = [];
    end
    totalEdepSum = scorerValues(fullfile(sourceDir, "TotalEDep_Target.csv"));
    if isempty(totalEdepSum)
        run.totalEdepMeVPerPrimary = NaN;
    else
        run.totalEdepMeVPerPrimary = totalEdepSum(1) ./ cfg.nPrimaries;
    end
    run.entranceFluenceSum = scorerValues(fullfile(sourceDir, "Fluence_TargetEntrance.csv"));
    run.exitFluenceSum = scorerValues(fullfile(sourceDir, "Fluence_TargetExit.csv"));
    run.gammaFluenceSum = scorerValues(fullfile(sourceDir, "GammaFluence_Exit.csv"));
    run.neutronFluenceSum = scorerValues(fullfile(sourceDir, "NeutronFluence_Exit.csv"));
end

function doseMap = doseMapFromCsv(fileName, cfg)
    values = scorerValues(fileName);
    if isempty(values)
        doseMap = [];
        return;
    end
    values = values(:) ./ cfg.nPrimaries;
    if numel(values) == 65
        doseMap = repmat(values.', 14, 1);
    elseif numel(values) == 14 * 65
        doseMap = reshape(values, [14, 65]);
    else
        doseMap = [];
    end
end

function [backForCombination, note] = orientBackDoseForCombination(frontDose, backDose)
    % TOPAS usually scores the back-beam run in fixed target coordinates, so
    % the back profile should already be high near the downstream end. If a
    % run instead looks like a second front profile, mirror it before adding.
    frontDose = frontDose(:);
    backDose = backDose(:);
    n = min(numel(frontDose), numel(backDose));
    frontDose = frontDose(1:n);
    backDose = backDose(1:n);

    directLikeFront = profileDistance(backDose, frontDose);
    mirrorLikeFront = profileDistance(flipud(backDose), frontDose);

    if mirrorLikeFront < directLikeFront
        backForCombination = backDose;
        note = "Back run already appears mirrored in target coordinates; no flip applied.";
    else
        backForCombination = flipud(backDose);
        note = "Back run resembled the front profile; flipped before front+back superposition.";
    end
end

function [combinedMap, note] = combineDoseMaps(front, back)
    if isempty(front.doseMapGyPerPrimary)
        combinedMap = [];
        note = "No 2D dose map was loaded.";
        return;
    end
    if isempty(back) || isempty(back.doseMapGyPerPrimary)
        combinedMap = front.doseMapGyPerPrimary + fliplr(front.doseMapGyPerPrimary);
        note = "Front map + mirrored-front estimate";
        return;
    end
    frontMap = front.doseMapGyPerPrimary;
    backMap = back.doseMapGyPerPrimary;
    if ~isequal(size(frontMap), size(backMap))
        combinedMap = frontMap + fliplr(frontMap);
        note = "Back map size mismatch; front map + mirrored-front estimate";
        return;
    end
    directLikeFront = profileDistance(backMap(:), frontMap(:));
    mirrorLikeFront = profileDistance(fliplr(backMap), frontMap);
    if mirrorLikeFront < directLikeFront
        backForCombination = backMap;
        note = "Front map + back map in target coordinates";
    else
        backForCombination = fliplr(backMap);
        note = "Front map + flipped back map";
    end
    combinedMap = frontMap + backForCombination;
end

function d = profileDistance(a, b)
    a = normalizeVector(a);
    b = normalizeVector(b);
    n = min(numel(a), numel(b));
    d = norm(a(1:n) - b(1:n));
end

function values = scorerValues(fileName)
    if ~isfile(fileName)
        values = [];
        return;
    end
    fid = fopen(fileName, "r");
    rawLines = textscan(fid, "%s", "Delimiter", "\n", "Whitespace", "");
    fclose(fid);
    rawLines = rawLines{1};
    rows = {};
    for i = 1:numel(rawLines)
        line = strtrim(rawLines{i});
        if isempty(line) || startsWith(string(line), "#")
            continue;
        end
        tokens = regexp(line, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
        nums = str2double(tokens);
        if ~isempty(nums)
            rows{end + 1} = nums(:).'; %#ok<AGROW>
        end
    end
    if isempty(rows)
        values = [];
        return;
    end
    widths = cellfun(@numel, rows);
    if isscalar(rows)
        values = rows{1}.';
    elseif all(widths >= 4)
        values = cellfun(@(r) r(end), rows).';
    else
        values = cellfun(@(r) r(1), rows).';
    end
end

function out = absoluteDose(doseGyPerPrimary, currentUA, mode, sides, cfg)
    ePerSecond = currentUA * 1.0e-6 / cfg.eChargeC;
    timeS = irradiationTimeSeconds(mode, sides, cfg);
    out = struct();
    out.doseRateGyS = doseGyPerPrimary .* ePerSecond;
    out.cumulativeGy = out.doseRateGyS .* timeS;
end

function t = irradiationTimeSeconds(mode, sides, cfg)
    if mode == "warm"
        hoursPerSide = cfg.warmHoursPerSide;
    else
        hoursPerSide = cfg.coldHoursPerSide;
    end
    t = hoursPerSide * 3600.0 * sides;
end

function fluence = deliveredFluence(currentUA, mode, sides, cfg)
    ePerSecond = currentUA * 1.0e-6 / cfg.eChargeC;
    fluence = ePerSecond * irradiationTimeSeconds(mode, sides, cfg) / cfg.interactionAreaCm2;
end

function heatW = heatLoadFromEDep(edepMeVPerPrimary, currentUA, cfg)
    ePerSecond = currentUA * 1.0e-6 / cfg.eChargeC;
    heatW = edepMeVPerPrimary * cfg.mevToJ * ePerSecond;
end

function nElectronsCm3 = electronDensityND3(densityGcm3, cfg)
    moleculesCm3 = densityGcm3 / cfg.nd3MolarMassGmol * cfg.avogadro;
    nElectronsCm3 = moleculesCm3 * cfg.nd3ElectronsPerMolecule;
end

function [rangeCm, rangeGcm2] = csdaRangeEstimate(eMeV, densityGcm3)
    a1 = 0.2335; a2 = 1.209; a3 = 1.078; a4 = 0.5842;
    rangeGcm2 = a1 * eMeV^a2 / (1.0 + a3 * eMeV * exp(-a4 * eMeV));
    rangeCm = rangeGcm2 / densityGcm3;
end

function s = yesNo(flag)
    if flag, s = "YES"; else, s = "NO"; end
end

function plotSpectrum(varargin)
    if nargin == 3
        ax = gca;
        eCenters = varargin{1};
        values = varargin{2};
        color = varargin{3};
    else
        ax = varargin{1};
        eCenters = varargin{2};
        values = varargin{3};
        color = varargin{4};
    end
    if isempty(values), return; end
    n = min(numel(eCenters), numel(values));
    plot(ax, eCenters(1:n), values(1:n), "LineWidth", 3.0, "Color", color);
end

function y = normalizeVector(x)
    x = x(:);
    finiteVals = x(isfinite(x));
    if isempty(finiteVals) || max(abs(finiteVals)) == 0
        y = zeros(size(x));
    else
        y = x ./ max(abs(finiteVals));
    end
end

function y = depthVector(values, referenceDepthVector)
    values = values(:);
    if isempty(values)
        y = [];
        return;
    end
    nRef = numel(referenceDepthVector);
    if numel(values) ~= nRef
        y = interp1(linspace(0, 1, numel(values)), values, ...
            linspace(0, 1, nRef), "linear", "extrap").';
    else
        y = values;
    end
end

function ax = nextPlotTab(tabGroup, tabTitle, pal)
    tab = uitab(tabGroup, "Title", tabTitle, "BackgroundColor", pal.page);
    ax = axes("Parent", tab, "Units", "normalized", "Position", [0.15 0.13 0.76 0.78]);
    ax.Color = pal.panel;
end

function saveAxesCopy(ax, outputDir, baseName)
    % Axes inside inactive uitabs are not always accepted by exportgraphics
    % on all MATLAB releases. Copy the axes into a hidden figure for export
    % so the user sees only the tabbed browser window.
    tempFig = figure("Visible", "off", "Color", ax.Parent.BackgroundColor, ...
        "InvertHardcopy", "off", "Position", [100 100 980 720]);
    tempAx = copyobj(ax, tempFig);
    set(tempAx, "Units", "normalized", "Position", [0.15 0.13 0.76 0.78]);
    oldColorbars = findall(ax.Parent, "Type", "ColorBar");
    if ~isempty(oldColorbars)
        for i = 1:numel(oldColorbars)
            try
                if oldColorbars(i).Axes == ax
                    tempCb = colorbar(tempAx);
                    tempCb.Label.String = oldColorbars(i).Label.String;
                    tempCb.Label.Color = oldColorbars(i).Label.Color;
                    tempCb.Color = oldColorbars(i).Color;
                end
            catch
            end
        end
    end
    exportgraphics(tempFig, fullfile(outputDir, baseName + ".png"), "Resolution", 300);
    close(tempFig);
end

function styleAxes(ax, pal)
    ax.Color = pal.panel;
    ax.XColor = pal.ink;
    ax.YColor = pal.ink;
    ax.GridColor = [0.60 0.68 0.82];
    ax.MinorGridColor = [0.82 0.86 0.94];
    ax.GridAlpha = 0.35;
    ax.LineWidth = 1.2;
    ax.FontSize = 12;
    ax.FontWeight = "bold";
    grid(ax, "on");
    box(ax, "on");
end

function styleLegend(lg, pal)
    lg.Color = [1.000 0.990 0.930];
    lg.TextColor = pal.ink;
    lg.EdgeColor = pal.reference;
    lg.LineWidth = 0.8;
end
