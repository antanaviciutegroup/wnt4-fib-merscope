Code for processing MERSCOPE transcript data with Baysor and Seurat.

## Overview

This repository contains a small analysis workflow for:

- reformatting detected transcript tables from MERSCOPE for Baysor
- running Baysor segmentation on transcript coordinates
- loading Baysor output into Seurat for clustering and marker analysis

## Contents

- `reformat_merscope_cell_ids_for_baysor.R`  
  Converts `detected_transcripts.csv` to `detected_transcripts_for_baysor.csv` by replacing unassigned `cell_id` values (`-1`) with `0`.

- `baysor.sh`  
  SLURM submission script for running Baysor on `detected_transcripts_for_baysor.csv`.

- `seurat_clusters_from_baysor.R`  
  Reads Baysor segmentation output, filters low-confidence assignments, builds a Seurat object, runs clustering/UMAP and assigns cluster labels.
