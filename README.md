# mqqc — MaxQuant-based LC-MS instrument QC

A scaffold for semi-automatic LC-MS instrument qualification in bottom-up proteomics,
built on standard MaxQuant output.

**Goal:** make the hand-off between instrument operators and data analysis efficient.
Operators keep working in the MaxQuant GUI they already know; this script turns each
new QC search into a consistent set of interactive plots and a summary table, and
accumulates results over time so instrument performance can be followed as a
**time course** rather than judged run by run.

## What it produces

Running `updateMQQC.R` in an instrument folder writes:

| Output | Content |
|---|---|
| `QC_IDs_lineplots.{pdf,html}` | Precursor / peptide / protein IDs per run over acquisition days |
| `QC_Mass_error_ppm.{pdf,html}` | Precursor mass error distributions per run |
| `QC_ID_density_along_RT.{pdf,html}` | Identification density along the gradient (latest *n* runs) |
| `QC_Retention_length.{pdf,html}` | Chromatographic peak widths |
| `QC_HeLa_3peptides_RTs.{pdf,html}` | Retention-time stability of three reference HeLa peptides |
| `QC_Summary_table.csv` | MaxQuant `summary.txt` metrics merged with run annotation |
| `QC_plots_history/` | Dated backups of the previous plots and annotation |

HTML versions are interactive (plotly) — hover a point to see the raw file.

## Folder layout

One copy of this repo per instrument; the folder name is used as the instrument name in
plot titles.

```
<InstrumentName>/
├── updateMQQC.R
├── mqpar.xml                 # MaxQuant parameter template (variants: mqpar_tims1.xml, mqpar_zeno1.xml)
├── runAnnotation.csv         # operator-maintained run metadata
├── mq_txt_results/
│   ├── txt_220815/           # one MaxQuant txt/ folder per QC batch: txt_<YYMMDD>
│   │   ├── evidence.txt
│   │   └── summary.txt
│   └── txt_<YYMMDD>/ ...
└── QC_plots_history/
```

## Usage

1. **Acquire** QC runs (e.g. HeLa digest) on the instrument.
2. **Search** them in the MaxQuant GUI, loading the matching `mqpar*.xml` as template.
3. **Copy** the resulting `combined/txt` folder to `mq_txt_results/txt_<YYMMDD>`.
4. **Run** the script once to generate `runAnnotation_template.csv`:
   ```bash
   Rscript updateMQQC.R
   ```
5. **Annotate**: copy new rows from `runAnnotation_template.csv` into
   `runAnnotation.csv` and replace the `FILL_IN` fields
   (operator, LC, column, LC method, sample, load, QC box).
6. **Run again** — plots and summary table are regenerated across all batches.

Steps 3–6 are the only manual work per QC batch.

## Requirements

R ≥ 4.0 with `data.table`, `ggplot2`, `plotly`, `ggrepel`, `htmlwidgets`:

```r
install.packages(c("data.table", "ggplot2", "plotly", "ggrepel", "htmlwidgets"))
```

MaxQuant (any recent version producing `evidence.txt` and `summary.txt`).

## Status

Working scaffold (parameter templates included for timsTOF and ZenoTOF); expect to adapt
paths, reference peptides and annotation fields to your lab.

## License

GPL-3.0 — see [LICENSE](LICENSE).

---
<sub>README optimized with Claude Opus 5.5.</sub>
