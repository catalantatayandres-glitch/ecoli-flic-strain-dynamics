# Results — Best vs Standard Pipeline Comparison

Direct sequence-level and biological comparison of the BEST pipeline
(`Strict · Binned · FALSE pooling`) against the STANDARD pipeline
(`Loose · PacBio · FALSE pooling`). Three complementary analyses assess
how different the two pipelines are and what the biological validity of
their respective ASV sets is.

---

## Figure 1 — Edit Distance Distributions

![Edit distance distribution](Figure1_Edit_Distance_Distribution.png)

Edit distance distributions to the closest matching sequence, shown for three
comparisons: within the BEST pipeline (blue, top), within the STANDARD pipeline
(green, middle), and cross-pipeline (orange, bottom). The dashed vertical line
separates the full-resolution region (0–5 bp) from the compressed region. The
x-axis uses a positional encoding so short distances occupy more visual space.

**Key findings:**
- The **within-BEST** distribution (panel 1) shows very few ASVs at edit distance 1, suggesting that DADA2 has effectively collapsed near-identical sequences under Strict filtering — few sequences that should have been merged were instead called as distinct variants
- The **within-STANDARD** distribution (panel 2) shows a substantially larger peak at edit distance 1 (~110 ASVs). In amplicon sequencing, ASVs differing by only 1 bp from another ASV in the same dataset are potentially denoising artefacts — sequences that should have been collapsed into their closest neighbour but were instead called as distinct variants due to insufficient error correction. This raises the possibility that a portion of the Standard pipeline's additional ASVs represent residual noise rather than true biological diversity, though further validation would be needed to confirm this
- The **cross-pipeline** distribution (panel 3) reveals a large peak at edit distance 0 (~260 ASVs), corresponding to sequences shared between both pipelines, followed by a secondary peak at distance 1. A substantial number of sequences have no identical counterpart in the other pipeline, indicating meaningful divergence in the ASV sets produced by each approach

---

## Figure 2 — Integrated Dendrogram

![Integrated dendrogram](Figure2_Integrated_Dendrogram.png)

Hierarchical clustering of all ASVs by sequence similarity (edit distance). Each row
is one ASV. Annotation tracks show, from left to right: pipeline origin (Best only =
purple, Shared = grey-blue, Standard only = yellow), biological group based on BLAST
trust tier assignment (E. coli = green, Ambiguous = orange, Non-E. coli = red),
scaled relative abundance (blue gradient heatmap), prevalence across samples (blue
gradient), and per-position nucleotide composition.

**Biological validation — trust tier system:**
Each ASV was submitted to NCBI BLAST and classified into one of five trust tiers:
- **Tier 1** — all hits show 100% identity and 100% query coverage to *E. coli*, or exact match to the WGS reference collection → highest confidence
- **Tier 2** — at least one perfect *E. coli* hit, but not all
- **Tier 3** — no perfect hits, but all hits relate to *E. coli*
- **Tier 4** — no perfect hits, mixed *E. coli* and non-*E. coli* taxonomy
- **Tier 5** — no *E. coli* hits at all

For downstream analyses, tiers were collapsed into three biological groups: **E. coli** (Tiers 1–2), **Ambiguous** (Tiers 3–4), and **Non-E. coli** (Tier 5).

**Key findings:**
- The most abundant and prevalent ASVs — visible as the darkest cells in the abundance and prevalence heatmaps — are predominantly classified as **E. coli** (green) and are **shared between both pipelines**
- ASVs unique to the Standard pipeline (yellow strip) are predominantly associated with low abundance, low prevalence, and enrichment for Non-*E. coli* and Ambiguous biological groups — consistent with the possibility that these represent noise or sequences of uncertain biological origin
- ASVs unique to the Best pipeline (purple strip) are very few in number and similarly show low abundance and prevalence
- The shared ASV pool has the highest proportion of biologically validated *E. coli* sequences, confirming that the core biological signal is robustly captured by both pipelines
- The nucleotide composition panel (right) reveals the sequence-level diversity underlying the clustering, with the dominant shared *E. coli* ASVs showing distinct, conserved per-position profiles

---

## Figure 3 — Bio Group Composition by ASV Origin

![Bio group by origin](Figure3_BioGroup_by_Origin.png)

Proportion (top) and absolute count (bottom) of biological groups (E. coli, Ambiguous,
Non-E. coli) for each ASV category: Unique to Standard, Shared, and Unique to Best.

**Key findings:**
- **Shared ASVs** show the highest proportion of biologically validated *E. coli* sequences (31.3%, n = 84 ASVs), confirming that the core signal captured by both pipelines is predominantly real *E. coli* biology
- **ASVs unique to the Standard pipeline** show a strikingly different composition: only 9.2% are classified as *E. coli* (n = 12), while 50.8% are Non-*E. coli* (n = 66) and 40% are Ambiguous (n = 52). Over 90% of Standard-exclusive ASVs therefore lack reliable biological validation as true *E. coli* sequences
- **ASVs unique to the Best pipeline** are negligible in absolute terms (n = 16 total), and are predominantly Non-*E. coli* and Ambiguous, consistent with the strict pipeline occasionally discarding low-confidence sequences
- These results reveal a clear **sensitivity–specificity trade-off** between the two pipelines: the Standard pipeline recovers a small amount of additional true *E. coli* diversity (12 extra validated ASVs) at the cost of generating substantially more sequences of uncertain or non-*E. coli* origin. The Best pipeline provides a cleaner and more conservative output, though the possibility that it discards a small number of genuinely rare *E. coli* variants cannot be excluded

---

*Figures generated by `04_BEST_vs_STANDARD.Rmd`. Full methods in the companion script and the project report.*
