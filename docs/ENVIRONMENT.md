# Environment and dependencies

Packaging and validation used R 4.4.3 on macOS arm64. Versions below describe the local validation environment; they are not a claim that every archived upstream analysis used the same versions. Run `Rscript scripts/check_environment.R` on another machine.

| Package | Installed version |
|---|---|
| BiocParallel | 1.40.2 |
| MCPcounter | 1.2.0 |
| Matrix | 1.7.4 |
| RColorBrewer | 1.1.3 |
| Seurat | 5.4.0 |
| copykat | 1.1.0 |
| cowplot | 1.2.0 |
| dplyr | 1.2.0 |
| fgsea | 1.32.4 |
| ggplot2 | 4.0.2 |
| ggpubr | 0.6.3 |
| gridExtra | 2.3 |
| gtable | 0.3.6 |
| harmony | 1.2.4 |
| immunedeconv | 2.1.0 |
| jsonlite | 2.0.0 |
| limma | 3.62.2 |
| msigdbr | 25.1.1 |
| patchwork | 1.3.2 |
| pheatmap | 1.0.13 |
| ragg | 1.5.0 |
| readr | 2.2.0 |
| readxl | 1.4.5 |
| reshape2 | 1.4.5 |
| survival | 3.8.6 |
| survminer | 0.5.1 |
| svglite | 2.2.2 |
| tibble | 3.3.1 |
| tidyr | 1.3.2 |
| writexl | 1.5.4 |

Base R packages `stats`, `utils`, `tools`, `grid` and `grDevices` are also used. The Git guard requires Python 3 (standard library only) and Git. Plot export requires Cairo support; Arial is the manuscript font, configurable via `font_family`.

Most packages are on CRAN. `limma`, `fgsea` and `BiocParallel` are Bioconductor packages; `MCPcounter`, `immunedeconv` and `copykat` use their respective project distributions. The upstream single-cell route also requires Harmony and sufficient RAM for the complete count matrix. Package installation is left explicit; the analysis scripts do not install packages or change the user's R library.

## External reference databases

The installed `msigdbr` uses its matching MSigDB release. The runner saves the gene-set data frame and database-version field outside Git; configure `hallmark_rds` to reuse the same snapshot. An R package version alone does not pin external database content.

MCPcounter's default human-symbol signature is fetched from the official project's public signature table at runtime. An internet connection is needed for that call in the installed package. The data are not copied into this code repository. For exact long-term reconstruction, archive the referenced signature version alongside the external inputs.

`fgsea` is explicitly given `BiocParallel::SerialParam()` so an analysis does not depend on machine-specific worker counts or local socket availability. CopyKAT's `workers` setting is independently configurable.

Full raw-count preprocessing/CopyKAT is optional and computationally expensive. Reuse the archived processed Seurat object and prediction files when reproducing the manuscript panels; use the optional upstream route to reconstruct inputs and separately review effects of software-version changes.
