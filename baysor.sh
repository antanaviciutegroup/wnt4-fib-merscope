#!/bin/sh
#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=20
#SBATCH --mem=200G
#SBATCH --partition=long

export JULIA_NUM_THREADS=20 

OUT=/PATH/TO/OUTPUT
IN=detected_transcripts_for_baysor.csv

baysor run $IN :cell_id -m 11 -s 5.55 -x global_x -y global_y -z global_z --save-polygons=geojson -p -o $OUT --n-clusters 12 --prior-segmentation-confidence 0.7 --count-matrix-format tsv