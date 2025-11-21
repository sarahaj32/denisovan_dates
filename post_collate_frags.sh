#!/bin/bash

#LASIDAD
sbatch ./collate_frags.sh /global/scratch/p2p3/pl1_moorjani/lauritsskov2/ArchaicSegments/Results/03decode/LASIDAD LASIDAD

#1000G
sbatch ./collate_frags.sh /global/scratch/p2p3/pl1_moorjani/lauritsskov2/ArchaicSegments/Results/03decode/1000genomes 1000G

#HGDP
sbatch ./collate_frags.sh /global/scratch/p2p3/pl1_moorjani/lauritsskov2/ArchaicSegments/Results/03decode/HGDP/HGDP HGDP

#GAv1 (I needed to rerun these to add on the -extrainfo)
# I match Yulin's results if I use the strict mask
sbatch ./collate_frags.sh data/GA100K/hmmix_strict/02decode/ GAv1
