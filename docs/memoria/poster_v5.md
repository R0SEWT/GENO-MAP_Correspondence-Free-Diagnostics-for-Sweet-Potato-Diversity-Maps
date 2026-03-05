# Poster A1 — GENO-MAP v5 (claim-aligned storyboard)

**Format:** A1 portrait (594 × 841 mm) · Frontiers style · 300 DPI

---

## Claims (the three things this poster proves)

| ID | Claim | Key evidence |
|----|-------|--------------|
| **A** | Correspondence-free validation is necessary | IDs disjuntos → no cross-panel alignment without manifest |
| **B** | PCA geometry is stable; neighborhoods degrade gracefully | SS ≥ 0.91 at 5% markers; $J_{\text{nbr}}$ monotonic; imputation insensitive |
| **C** | AE is not operationally justified in $n \ll p$ | Trust gains marginal; stability gap ~0.46; crossover $n^*$ irrelevant |

---

## Figure set (final — claim-mapped)

| Fig | Content | Claim | Role | Size | Status |
|-----|---------|-------|------|------|--------|
| 1 | Disjoint IDs diagram | A | Claim-bearing | Small | 🟡 MOCK (draw.io) |
| 2 | Pipeline flowchart | — | Supporting | Small | 🟡 MOCK (draw.io) |
| 3 | Panel geometry + QA thresholds | B | Claim-bearing | Medium | ✅ `fig3_panel_geometry.png` |
| 4 | Robustness HERO (sub + miss + SS) | B | **Claim-bearing (HERO)** | **XL** | ✅ `fig4_robustness_hero.png` |
| 5* | UMAP scatter (1 panel only) | — | Decorative/inset | Mini | ✅ `fig5_umap_single.png` |
| 6 | PCA vs AE trade-off bars | C | Claim-bearing | Medium | ✅ `fig6_pca_vs_ae.png` |

**Dropped from poster** (→ supplementary via QR):
- Fig 5 full (4-panel UMAP) — decorative, doesn't support claims
- Fig 10 (stability frontier) — secondary evidence for C, "fine detail"

---

## Layout (A1 portrait)

```
┌──────────────────────────────────────────────────────────┐
│                                                          │
│  GENO-MAP: Correspondence-Free Diagnostics and           │
│  Robustness Curves for Sweet Potato Diversity Maps       │
│                                                          │
│  Rody Vilchez · UNALM / CIP                              │
│                                                          │
│  ┌────────────────────────────────────────────────────┐  │
│  │ Validating ultra-wide DArT panels without shared   │  │
│  │ IDs: PCA geometry is robust, neighborhoods degrade │  │
│  │ gracefully, and autoencoders add complexity without │  │
│  │ stability in the n≪p regime.                       │  │
│  └────────────────────────────────────────────────────┘  │
│                                                          │
├──────────────────────────┬───────────────────────────────┤
│                          │                               │
│  1. WHY CORRESPONDENCE-  │  2. METHOD                    │
│     FREE?                │                               │
│                          │  [FIG 2: pipeline]            │
│  [FIG 1: IDs disjuntos]  │   DArT → impute → PCA →      │
│                          │   UMAP → kNN → JSON           │
│  • Global: CIP accession │                               │
│    codes (n = 5 970)     │  • PCA-30D = analytic space   │
│  • LowDensity: plate/    │  • UMAP = visualization only  │
│    well DArT (n = 630)   │  • kNN (k=15, cosine)         │
│  • No manifest → no      │  • Validation: geometry diag  │
│    alignment possible    │    + robustness curves         │
│                          │                               │
│  [FIG 5*: UMAP inset]    │  Mini-table:                  │
│   Global SNP example     │  ┌────────┬────┬─────┬──────┐ │
│                          │  │Panel   │ n  │  p  │ n/p  │ │
│                          │  ├────────┼────┼─────┼──────┤ │
│                          │  │G-SNP   │5970│ 20k │0.297 │ │
│                          │  │G-Silico│5970│ 57k │0.103 │ │
│                          │  │LD-SNP  │ 630│ 62k │0.010 │ │
│                          │  │LD-Sil. │ 635│ 38k │0.017 │ │
│                          │  └────────┴────┴─────┴──────┘ │
│                          │                               │
├──────────────────────────┴───────────────────────────────┤
│                                                          │
│  3. RESULTS — ROBUSTNESS (HERO)                          │
│                                                          │
│  ┌────────────────────────────────────────────────────┐  │
│  │                                                    │  │
│  │  [FIG 4: Robustness HERO — ~35% of poster area]    │  │
│  │                                                    │  │
│  │  (A) Marker subsampling: J_nbr from 0.43 to 0.90  │  │
│  │      Inset: SS ≥ 0.91 even at 5% markers          │  │
│  │  (B) Missing injection: monotonic degradation      │  │
│  │      J_nbr = 0.77–0.86 at +20% MCAR               │  │
│  │                                                    │  │
│  └────────────────────────────────────────────────────┘  │
│                                                          │
├──────────────────────────┬───────────────────────────────┤
│                          │                               │
│  4. PANEL GEOMETRY       │  5. PCA vs AUTOENCODER        │
│                          │                               │
│  [FIG 3: dot-plot with   │  [FIG 6: grouped bars         │
│   QA threshold lines]    │   trust + stability]          │
│                          │                               │
│  • No 1D collapse        │  • Trust: AE +1 pp (Global),  │
│    (PC1 < 11%)           │    −8 pp (LowDensity)         │
│  • Eff. rank 10–16       │  • Stability: PCA 0.88–0.91   │
│  • Reciprocity 0.65–0.70 │    vs AE 0.16–0.52            │
│  • Flags: EXTREME-WIDE,  │  • "Marginal trust does not   │
│    DISCONNECTED, HIGH-   │    justify instability"        │
│    MISSINGNESS           │                               │
│                          │                               │
├──────────────────────────┴───────────────────────────────┤
│                                                          │
│  KEY TAKEAWAYS                                           │
│                                                          │
│  1. PCA subspace geometry is robust: SS ≥ 0.91 under     │
│     95% marker dropout — structure is genome-wide,       │
│     not marker-specific.                                 │
│                                                          │
│  2. kNN neighborhoods degrade monotonically and          │
│     predictably under perturbation — no cliff effects.   │
│                                                          │
│  3. Correspondence-free validation (geometry + robustness│
│     curves) provides rigorous QA without shared IDs.     │
│                                                          │
│  4. Autoencoders offer ≤1 pp trust gain but 2–5×         │
│     worse topological stability — not justified for      │
│     operational use in n≪p sweet potato panels.          │
│                                                          │
│  ┌─────────┐  GitHub: [QR]  CIP Dataverse: [QR]         │
│  │  [QR]   │  Contact: rvilchez@lamolina.edu.pe          │
│  └─────────┘                                             │
│                                                          │
└──────────────────────────────────────────────────────────┘
```

---

## One-liner claim (below title)

> **Validating ultra-wide DArT panels without shared IDs: PCA geometry is robust, neighborhoods degrade gracefully, and autoencoders add complexity without stability in the $n \ll p$ regime.**

---

## Key takeaways (bottom of poster, 4 bullets)

1. **PCA subspace geometry is robust**: subspace similarity ≥ 0.91 under 95% marker dropout — structure is genome-wide, not marker-specific.

2. **kNN neighborhoods degrade monotonically** and predictably under perturbation (marker subsampling, MCAR injection) — no cliff effects, enabling quality thresholds.

3. **Correspondence-free validation** (geometry diagnostics + robustness curves) provides rigorous panel QA without shared identifiers — directly actionable for germplasm curators.

4. **Autoencoders offer ≤ 1 pp trustworthiness gain** but 2–5× worse topological stability — not operationally justified in the $n \ll p$ regimes observed in sweet potato DArT panels.

---

## Text budget per section (A1 constraint)

| Section | Max words | Content |
|---------|-----------|---------|
| Title + one-liner | 35 | As above |
| Why corr-free (left col) | 40 | 3 bullets: Global IDs, LD IDs, no manifest |
| Method (right col) | 40 | Pipeline steps + PCA=analytic, UMAP=viz |
| Mini-table | — | 4 rows, 4 cols (Panel, n, p, n/p) |
| Robustness hero caption | 50 | Fig 4 caption: sub + miss + SS inset |
| Panel geometry caption | 30 | Fig 3 caption: 3 metrics + QA flags |
| PCA vs AE caption | 40 | Fig 6 caption: trust + stability bars |
| Key takeaways | 80 | 4 bullets as above |
| **Total** | **~315** | Target ≤ 350 words for A1 |

---

## What was killed (lives in paper/repo, not on poster)

- Tabla 2–10 completas → paper
- Arquitectura AE detallada (ResBlock, dropout, noise, epochs) → paper §3.6
- Discusión larga (§5.1–5.6) → paper §5
- Limitaciones extensas (§7) → paper §7
- Trabajo futuro (§8) → paper §8
- Fig 10 (stability frontier) → supplementary PDF via QR
- Fig 5 full (4-panel UMAP) → supplementary PDF via QR
- Imputation comparison table → paper §4.5.3

---

## Supplementary figures (accessible via QR)

These figures are generated by the same script but NOT printed on the A1 poster:

| Fig | Content | File |
|-----|---------|------|
| 5 (full) | UMAP scatter, 4 panels | `fig5_umap_scatter.png` |
| 10 | Stability frontier (crossover $n^*$ + gap) | `fig10_stability_frontier.png` |

Both support Claim C but are "fine detail" that competes with the 10-second reading rhythm of a poster.

---

## Assets checklist

### Generated (in `docs/figures/poster/`)

- [x] `fig3_panel_geometry.{png,pdf}` — Claim B
- [x] `fig4_robustness_hero.{png,pdf}` — Claim B (HERO)
- [x] `fig5_umap_single.{png,pdf}` — decorative inset (1 panel)
- [x] `fig5_umap_scatter.{png,pdf}` — supplementary (4 panels)
- [x] `fig6_pca_vs_ae.{png,pdf}` — Claim C
- [x] `fig10_stability_frontier.{png,pdf}` — supplementary (Claim C detail)

### Manual (draw.io / Inkscape / TikZ)

- [ ] Fig 1 — Disjoint IDs conceptual diagram (Claim A)
- [ ] Fig 2 — Pipeline flowchart (supporting)
- [ ] QR codes (GitHub repo + CIP Dataverse)

---

*Compiled: March 2026*
*Script: `scripts/generate_poster_figures.py`*
*Output: `docs/figures/poster/`*
