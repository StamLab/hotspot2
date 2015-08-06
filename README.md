# README #

BEDOPS enabled Hotspot.
The software is a fast implementation of the hotspot algorithm for identifying regions of enrichment
in next generation sequencing assays.

### Requirements ###

* BEDOPS >2.4 (http://bedops.readthedocs.org/)
* Python 2.7 (w/ numpy)
* random-lines (from https://bitbucket.org/jvierstra/bio-tools)

### Running ###

The software is modular -- various tasks will utlize somewhat different scripts.

Example on a single core:

	# simulate data for FDR calculation
	bash /home/jvierstra/proj/code/hotspot/scripts/simulate.cutcounts.sh --tmpdir=/tmp --starch-output fragments.starch ${GENOME_MAPPABILITY_FILE} cutcounts.simulated.starch

	# Call hotspots on observed data
	bash ${HOTSPOT_DIR}/scripts/hotspot.run.sh --tmpdir=/tmp cutcounts.starch ${GENOME_MAPPABILITY_FILE} hotspots.unthresholded.bed

	# Call hotspots on simulated data
	bash ${HOTSPOT_DIR}/scripts/hotspot.run.sh --tmpdir=/tmp cutcounts.simulated.starch ${GENOME_MAPPABILITY_FILE} hotspots.simulated.bed

	# Compute Q-values 
	python ${HOTSPOT_DIR}/scripts/compute_q_values.py hotspots.simulated.bed hotspots.unthresholded.bed > hotspots.fdr.bed

	# Select threshold and merge
	bash ${HOTSPOT_DIR}/scripts/hotspot.thresh-merge.sh --q-thresh=0.01 hotspots.fdr.bed hotspots.fdr.0.01.bed 

Example to run chromosomes on separate processes (e.g., cluster):

	for i in $(seq 1 22); do
		chrom="chr$i"
		bash ${HOTSPOT_DIR}/scripts/hotspot.run.sh --tmpdir=/tmp --contig=$chrom cutcounts.starch ${GENOME_MAPPABILITY_FILE} hotspots.${chrom}.unthresholded.bed
	done

* The "--contig=" option works on both the hotspot and simulate tags scripts. 

### TODO ###

* Badspot implementation
* Change downsample script to downsample from fragments, not cutcounts
