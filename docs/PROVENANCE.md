# Code provenance

Sources were read and adapted on 2026-09-30; originals remain unchanged. Paths are relative to the original research workspace, not runtime dependencies.

| Original source | SHA-256 |
|---|---|
| `20260721_Figures_redraw/Figure1.R` | `5457e5c19324c7274d4b068b5bb56e3dbbf84659af9a4d4166db6e5d791b02eb` |
| `20260721_Figures_redraw/Figure2.R` | `0f697dc25638e5027c47a0ab822fa9c9f3076e55bff8e18e9bb48fe54563e19c` |
| `Data/BRCA_scRNASeq/20260311_New_ScRNASeq_Results/01_QC_and_Annotation.R` | `d5049bda3ba42c25ccf554d2c8f80fcdb318d1159603177bfe3130905e9b4e73` |
| `Data/BRCA_scRNASeq/20260311_New_ScRNASeq_Results/03_OSM_Response_Signature.R` | `e545818ea7f41efcb9aa2493ce82bea0027f4b60c9903a48e19239f18577d411` |
| `20260111_New_Experiments/01_preprocessing_qc.R` | `29215cd767911def6340d28cf9cce136db069ee6a1bc02929d4d7477a6884120` |
| `20260111_New_Experiments/11_preprocessing_qc_OE.R` | `356827b143c427d44452930c7ae85b3b52e91ae21ca30b94094aaa47a0be35e5` |
| `20260111_New_Experiments/02_differential_expression.R` | `41f6e279ae551177f4fd912e345f85478b281c2c62c6dee8d22df5d514df1a76` |
| `20260111_New_Experiments/12_de_analysis_OE.R` | `62a491bd4ff92f3ea6f8367cdcb486253766d8853504190aedf7e7c8dfc74912` |
| `Publication/20260911_JECCR_Manuscript/20260923_Figure4_Bridge/01_bridge_analysis.R` | `773b9690ecdc9c548bde37ad082a8b68c0228860c2ea7a1f51e3de19bb7c06be` |
| `Publication/20260911_JECCR_Manuscript/20260930_figure4e_spacing/build_figure4.R` | `40d024e8bcc9bbddf6c76533d42c111637492f928452277718bbee70f9c27a4e` |
| `Publication/20260921_Figure4_Baseline_Context_Mock/v4_baseline_stats_response_genes/01_build_mock_v4.R` | `97dbc80a116912069bc428aec475f2146c563105e998ec2a1b49456e72bc2236` |

## Adaptation boundaries

Absolute workspace paths were replaced with configuration. Final panel letters replace historical letters. Only analyses required for the specified PDF panels and their statistical correction families were retained. The Figure 4 bandwidth dependency was replaced with the original calculation from TPMs, and figure data are regenerated outside this repository.

Figure 1d uses MCP-counter explicitly; unrelated method fallbacks are removed. Its original local/global color-order assignment was repaired so the intended plot runs. Hallmark rank scores keep their original formula and are named as such rather than as GSVA/ssGSEA.

The user confirmed the corrected `20260721_Figures_redraw/Figure2.R` master as the authoritative source for Figure 1f–j. Its `run_fig2c()` supplies the Figure 1h calculation, with only function/output names adapted to final panel numbering and a computed-score CSV export added. The earlier all-lineage CopyKAT implementation is superseded and is not included as an alternative execution path. Panels h–j consistently retain the master’s epithelial-only label transfer.
