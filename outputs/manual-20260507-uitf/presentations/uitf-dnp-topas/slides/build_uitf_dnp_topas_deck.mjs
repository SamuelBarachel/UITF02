import { Presentation, PresentationFile } from "/Users/takwirira/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";
import { readFile, writeFile } from "node:fs/promises";

const root = "/Users/takwirira/Desktop/TOPAS/UITF";
const assetDir = `${root}/figures_matlab/pptx_deck_assets`;
const outDir = `${root}/outputs/manual-20260507-uitf/presentations/uitf-dnp-topas/output`;
const previewDir = `${root}/outputs/manual-20260507-uitf/presentations/uitf-dnp-topas/preview`;
const finalPptx = `${outDir}/UITF_DNP_TOPAS_MC_Report.pptx`;

const W = 1280;
const H = 720;
const C = {
  bg: "#F4F8FC",
  ink: "#071333",
  muted: "#516176",
  line: "#D8E2EE",
  blue: "#0B72D9",
  green: "#009D71",
  orange: "#E85D04",
  purple: "#7A2EF2",
  paleBlue: "#EAF4FF",
  paleGreen: "#E9F8F0",
  paleOrange: "#FFF0E3",
  white: "#FFFFFF",
  dark: "#071333",
};

const img = {
  super: `${assetDir}/01_front_back_superimposed_cold_dose.png`,
  dose: `${assetDir}/02_combined_cold_dose_depth.png`,
  map: `${assetDir}/03_2d_dose_map_superimposed.png`,
  fluence: `${assetDir}/04_electron_fluence_depth_unnormalized.png`,
  density: `${assetDir}/05_material_electron_density_unnormalized.png`,
  uniform: `${assetDir}/06_lateral_dose_uniformity_raw.png`,
  spectra: `${assetDir}/07_energy_spectra.png`,
  heat: `${assetDir}/08_heat_load_vs_current.png`,
  acronyms: `${assetDir}/09_acronyms_and_interpretation_notes.png`,
};

const presentation = Presentation.create();

function slide() {
  const s = presentation.slides.add();
  s.setViewportSize(W, H);
  s.background.fill = { type: "solid", color: C.bg };
  return s;
}

function box(s, x, y, w, h, fill = C.white, line = C.line) {
  return s.shapes.add({
    geometry: "roundRect",
    position: { left: x, top: y, width: w, height: h },
    fill: { type: "solid", color: fill },
    line: { style: "solid", fill: line, width: 1 },
  });
}

function textBox(s, text, x, y, w, h, opts = {}) {
  const sh = s.shapes.add({
    geometry: "rect",
    position: { left: x, top: y, width: w, height: h },
    fill: { type: "solid", color: opts.fill ?? C.bg },
    line: { style: "solid", fill: opts.line ?? (opts.fill ?? C.bg), width: opts.lineWidth ?? 0 },
  });
  sh.text = text;
  sh.text.typeface = opts.font ?? "Aptos";
  sh.text.fontSize = opts.size ?? 28;
  sh.text.color = opts.color ?? C.ink;
  sh.text.bold = opts.bold ?? false;
  sh.text.alignment = opts.align ?? "left";
  sh.text.verticalAlignment = opts.valign ?? "top";
  sh.text.wrap = "square";
  return sh;
}

function title(s, eyebrow, headline, sub = "") {
  textBox(s, eyebrow, 52, 30, 500, 24, { size: 16, bold: true, color: C.orange });
  textBox(s, headline, 52, 58, 930, 88, { size: 36, bold: true, color: C.ink });
  if (sub) textBox(s, sub, 52, 138, 950, 40, { size: 16, color: C.muted });
  s.shapes.add({
    geometry: "rect",
    position: { left: 52, top: 190, width: 1160, height: 2 },
    fill: { type: "solid", color: C.line },
    line: { style: "solid", fill: C.line, width: 0 },
  });
}

function footer(s, n) {
  textBox(s, "TOPAS / Geant4 MC | UITF DNP target irradiation | generated May 7, 2026", 52, 684, 760, 18, {
    size: 11,
    color: C.muted,
  });
  textBox(s, String(n).padStart(2, "0"), 1160, 678, 70, 28, {
    size: 14,
    bold: true,
    color: C.muted,
    align: "right",
  });
}

function metric(s, label, value, x, y, w, color = C.blue) {
  box(s, x, y, w, 86, C.white, C.line);
  textBox(s, label, x + 18, y + 14, w - 36, 20, { size: 13, bold: true, color: C.muted, fill: C.white });
  textBox(s, value, x + 18, y + 38, w - 36, 34, { size: 24, bold: true, color, fill: C.white });
}

function bullets(s, items, x, y, w, h, opts = {}) {
  const text = items.map((d) => `- ${d}`).join("\n");
  textBox(s, text, x, y, w, h, { size: opts.size ?? 18, color: opts.color ?? C.ink, fill: opts.fill ?? C.bg });
}

function imageCard(s, path, x, y, w, h, label = "") {
  box(s, x, y, w, h, C.white, C.line);
  if (label) textBox(s, label, x + 18, y + 14, w - 36, 22, { size: 14, bold: true, color: C.muted, fill: C.white });
  const topPad = label ? 44 : 18;
  s.images.add({
    path,
    alt: label || "Simulation figure",
    position: { left: x + 18, top: y + topPad, width: w - 36, height: h - topPad - 18 },
    fit: "contain",
  });
}

// 1. Title / thesis
{
  const s = slide();
  textBox(s, "UITF DNP Target Irradiation", 52, 42, 620, 26, { size: 18, bold: true, color: C.orange });
  textBox(s, "TOPAS Monte Carlo supports the proposal front/back irradiation protocol", 52, 78, 840, 104, {
    size: 38,
    bold: true,
  });
  textBox(
    s,
    "8 MeV electrons through the cryostat, aluminum spreading window, liquid helium bath, and packed ND3 target volume. V2 production is defined as 10M front histories plus 10M back histories.",
    52,
    198,
    680,
    70,
    { size: 17, color: C.muted },
  );
  metric(s, "Histories", "10M + 10M", 52, 318, 245, C.purple);
  metric(s, "Cold current", "1 uA", 322, 318, 190, C.blue);
  metric(s, "Warm current", "10 uA", 537, 318, 190, C.orange);
  metric(s, "Target", "250 mg ND3", 752, 318, 250, C.green);
  imageCard(s, img.super, 724, 70, 496, 220, "front/back dose proof");
  box(s, 52, 438, 1168, 160, C.paleBlue, "#C9DDF4");
  textBox(s, "Main claim", 82, 464, 160, 24, { size: 15, bold: true, color: C.blue, fill: C.paleBlue });
  textBox(
    s,
    "The simulation now matches the proposal's required two-sided irradiation logic: one pass is range-limited, while front + back superposition gives the target-volume dose field used for interpretation.",
    82,
    498,
    1040,
    52,
    { size: 24, bold: true, color: C.ink, fill: C.paleBlue },
  );
  footer(s, 1);
}

// 2. Proposal protocol and MC setup
{
  const s = slide();
  title(s, "PROPOSAL RECREATION", "The model tracks the beamline and batch protocol the proposal actually describes");
  const rows = [
    ["Beam", "8 MeV electrons, 1-10 uA"],
    ["Target", "250 mg packed ND3 in 6.5 cm x 1.42 cm volume"],
    ["Geometry", "cryostat windows, liquid helium bath, target holder, front/back source"],
    ["Warm run", "1.0e17 e-/cm2 per side at 10 uA, 2.85 h/side"],
    ["Cold run", "5.0e15 e-/cm2 per side at 1 uA, 1.42 h/side"],
    ["Protocol", "two-sided irradiation because 8 MeV electrons are range-limited"],
  ];
  rows.forEach((r, i) => {
    const y = 226 + i * 58;
    box(s, 60, y, 460, 44, i % 2 ? C.white : C.paleBlue, C.line);
    textBox(s, r[0], 78, y + 12, 130, 18, { size: 14, bold: true, color: C.blue, fill: i % 2 ? C.white : C.paleBlue });
    textBox(s, r[1], 205, y + 10, 295, 22, { size: 14, color: C.ink, fill: i % 2 ? C.white : C.paleBlue });
  });
  imageCard(s, img.super, 560, 226, 620, 306, "front/back run pairing");
  box(s, 560, 554, 620, 78, C.paleOrange, "#F6D5B6");
  textBox(
    s,
    "Single-bead microdosimetry was removed from the report figures: the proposal mentions 1 mm target spheres, but it does not use single-bead dose as an experimental observable.",
    584,
    576,
    560,
    38,
    { size: 16, color: C.ink, fill: C.paleOrange },
  );
  footer(s, 2);
}

// 3. Dose superposition
{
  const s = slide();
  title(s, "DOSE RESULT", "Front + back superposition is the key physics result");
  imageCard(s, img.super, 58, 220, 575, 390, "individual sides plus combined profile");
  imageCard(s, img.dose, 660, 220, 560, 390, "combined cold dose-depth profile");
  box(s, 70, 604, 1120, 48, C.paleGreen, "#C9EAD9");
  textBox(
    s,
    "Cold mean cumulative dose: one side 1.317e6 Gy; front + back 2.633e6 Gy. The combined curve is not just a copy of one run.",
    94,
    618,
    1040,
    22,
    { size: 18, bold: true, color: C.ink, fill: C.paleGreen },
  );
  footer(s, 3);
}

// 4. Spatial coverage
{
  const s = slide();
  title(s, "SPATIAL COVERAGE", "The combined dose map is now superimposed across depth and radius");
  imageCard(s, img.map, 54, 218, 705, 410, "2D R-Z dose map, front + back");
  imageCard(s, img.uniform, 790, 218, 380, 410, "raw radial dose profile");
  textBox(
    s,
    "This is the slide to use when someone asks whether the target volume is being covered, not just whether the central depth profile looks good.",
    80,
    632,
    1040,
    26,
    { size: 16, bold: true, color: C.muted },
  );
  footer(s, 4);
}

// 5. Electron quantities
{
  const s = slide();
  title(s, "ELECTRON TRANSPORT", "Fluence and material electron density are shown without normalization");
  imageCard(s, img.fluence, 58, 216, 560, 390, "electron fluence scorer");
  imageCard(s, img.density, 654, 216, 560, 390, "ND3 material electron density");
  box(s, 70, 612, 1110, 52, C.paleBlue, "#C9DDF4");
  textBox(
    s,
    "The density curve is not particle fluence. It is the ND3 electron density computed from composition and density; the packed target volume is much lower than solid ND3 because the 250 mg batch is distributed over the holder volume.",
    92,
    626,
    1060,
    26,
    { size: 15, color: C.ink, fill: C.paleBlue },
  );
  footer(s, 5);
}

// 6. Beam quality and thermal feasibility
{
  const s = slide();
  title(s, "FEASIBILITY CHECKS", "Spectra and heat load show whether the proposed currents are plausible");
  imageCard(s, img.spectra, 58, 216, 560, 390, "entrance electrons and exit photons");
  imageCard(s, img.heat, 654, 216, 560, 390, "target heat load versus current");
  bullets(
    s,
    [
      "Cold irradiation at 1 uA is the relevant low-heat condition for 2 K running.",
      "Warm 10 uA running is evaluated separately because heat load is not simultaneous front + back.",
      "Spectra confirm the nominal 8 MeV entrance beam and the secondary photon component.",
    ],
    90,
    620,
    1040,
    54,
    { size: 14, color: C.muted },
  );
  footer(s, 6);
}

// 7. What is presentation-ready
{
  const s = slide();
  title(s, "INTERPRETATION", "What this deck can safely claim");
  const claims = [
    ["Strong", "Geometry, beam energy, target dimensions, front/back timing, and fluence targets match the proposal."],
    ["Strong", "The v2 10M front/back analysis is loaded from both completed output folders and superimposed in depth and R-Z maps."],
    ["Useful", "Electron fluence and material electron density are unnormalized and separated conceptually."],
    ["Useful", "Heat load is quoted per irradiation side, matching the actual procedure."],
    ["Not claimed", "No measured DNP polarization or ESR radical population is predicted by TOPAS alone."],
  ];
  claims.forEach((c, i) => {
    const y = 226 + i * 72;
    const fill = i < 2 ? C.paleGreen : i < 4 ? C.paleBlue : C.paleOrange;
    const accent = i < 2 ? C.green : i < 4 ? C.blue : C.orange;
    box(s, 86, y, 1060, 54, fill, C.line);
    textBox(s, c[0], 110, y + 15, 120, 20, { size: 15, bold: true, color: accent, fill });
    textBox(s, c[1], 240, y + 12, 860, 24, { size: 16, color: C.ink, fill });
  });
  footer(s, 7);
}

// 8. Acronym appendix
{
  const s = slide();
  title(s, "APPENDIX", "Acronym quick reference for the report");
  imageCard(s, img.acronyms, 120, 214, 1040, 420, "definitions used in the figures");
  footer(s, 8);
}

const pendingImages = presentation.getPendingImageHydrationRequests();
const hydratedImages = await Promise.all(
  pendingImages.map(async (request) => ({
    assetId: request.assetId,
    contentType: request.contentType || "image/png",
    data: await readFile(request.uri),
  })),
);
presentation.hydrateImageAssets(hydratedImages);

const pptx = await PresentationFile.exportPptx(presentation);
await pptx.save(finalPptx);

for (let i = 0; i < presentation.slides.count; i += 1) {
  const s = presentation.slides.getItem(i);
  const png = await s.export({ format: "png", scale: 1 });
  await writeFile(`${previewDir}/slide_${String(i + 1).padStart(2, "0")}.png`, Buffer.from(await png.arrayBuffer()));
}

console.log(finalPptx);
