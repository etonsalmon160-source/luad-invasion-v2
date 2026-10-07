#!/usr/bin/env python3
# fig_compose_main.py —— 拼正文大图（版式见 compose.py）
from compose import compose

compose(rows=[[("P16_cohort_design", 1.0)],
              [("P1c_UMAP_stage_split", 1.0)]],
        fname="Fig1_cohort_and_design", width=7.0)

compose(rows=[[("P1a_UMAP_lineage", 1.0), ("P2a_lineage_composition_by_stage", 1.0)],
              [("P2c_UMAP_L2_subtype", 1.0), ("P2b_L2_composition_by_stage", 1.0)]],
        fname="Fig2_singlecell_atlas", width=7.0, wspace=0.10)

compose(rows=[[("P3b_RCTD6_composition_by_stage", 1.0), ("P3a_RCTD39_composition_by_stage", 1.0)]],
        fname="Fig3_spatial_deconvolution", width=7.0, wspace=0.18)

compose(rows=[[("P5_spatial_matrix", 1.0)]],
        fname="Fig4_spatial_matrix", width=7.4)

compose(rows=[[("P4a_spatial_domain_archetype", 1.0)],
              [("P6a_domain_RCTD_lineage_composition", 1.0), ("P7_domain_marker_heatmap", 0.62)],
              [("P8_domain_dotplot", 1.0)]],
        fname="Fig5_domain_identity", width=7.4, wspace=0.08)

compose(rows=[[("P15_depth_density_coupling", 1.0), ("P14_domain_spatial_cnv", 1.0)]],
        fname="Fig6_depth_density_coupling", width=7.4, wspace=0.12)

compose(rows=[[("P9a_domain_forest_continuous", 1.0)],
              [("P9b_domain_forest_lowhigh", 1.0)]],
        fname="Fig7_domain_prognosis", width=7.0)

compose(rows=[[("P10_scmg_workflow", 1.0)],
              [("P11_scmg_positive_control", 0.62), ("P20_signature_qc", 1.0)]],
        fname="Fig8_reversal_pipeline", width=7.4, wspace=0.10)

compose(rows=[[("P12_compound_reversal_singlecell", 1.0)],
              [("P18_method_agreement", 1.0)],
              [("P17_split_half_stability", 1.0)],
              [("P13_compound_reversal_crossmodal", 1.0)],
              [("P19_domain_reversal", 1.0)]],
        fname="Fig9_reversal_candidates", width=7.4)
