###########################################################
# You may be able to get this to work with shell          #
# I tend to lean bash and have used that environment here #
###########################################################
SHELL := /bin/bash

###########################
# keep intermediate files #
###########################
.SECONDARY:

#####################################
# Key macros used extensively below #
# and values loaded from files      #
#####################################
RELEASE := 95
PY_CPU_FILE := changeable_parameters/num_motif_cpus.txt
MULTI_CPUS := $(shell cat ${PY_CPU_FILE})
SPECIES_FILE := changeable_parameters/species.txt
SPECIES := $(shell cat ${SPECIES_FILE})

##############
# CasX motif #
##############
SEARCH_MOTIF := TTCNNNNNNNNNNNNNNNNNNNNN

#################################################################
# Try to keep directory clean with subdirs for specific outputs #
#################################################################
DIRS := casx_sites_out gc_content logs gtfs fasta R_analysis figures
GTFFILES := $(addprefix gtfs/, $(addsuffix .gtf, $(SPECIES)))
FASTAFILES := $(addprefix fasta/, $(addsuffix .fa, $(SPECIES)))
GCFILES := $(addprefix gc_content/, $(addsuffix _gcPerc.txt, $(SPECIES)))
CASX_SITES_OUT := $(addprefix casx_sites_out/, $(addsuffix _casx_sites.csv.gz, $(SPECIES)))
ANNOTATED_SITES := $(addprefix R_analysis/, $(addsuffix _annotated_Casx_sites.csv, $(SPECIES)))
CUTSPERGENE := $(addprefix R_analysis/, $(addsuffix _cutsPerGene.csv, $(SPECIES)))
ROUTPUT := R_analysis/Casx_PAM_frequency.csv.gz R_analysis/Casx_site_counts.csv.gz $(ANNOTATED_SITES) $(CUTSPERGENE)
MOUSE_TARGETS := R_analysis/mouse_target_upset.csv.gz
HAMMING := R_analysis/mouse_hamming_distance.csv

##########################
# main Make run function #
##########################
all: $(DIRS) $(GCFILES) $(CASX_SITES_OUT) $(ROUTPUT) $(MOUSE_TARGETS) $(HAMMING)
	@echo Macros used in this makefile
	@echo DIRS: $(DIRS)
	@echo FASTA: $(FASTAFILES)
	@echo GC: $(GCFILES)
	@echo GTF: $(GTFFILES)
	@echo SPECIES: $(SPECIES)
	@echo motif scan multicpu: $(MULTI_CPUS)
	@echo CasX site files: $(CASX_SITES_OUT)
	@echo R output: $(ROUTPUT)

###############
# Directories #
###############
$(DIRS):
	mkdir -p $@

#######################################
# Calculate GC content of each genome #
#######################################
gc_content/%_gcPerc.txt: fasta/%.fa
	python ./fasta_gc_content_calc.py $< --output_file $@ 1>logs/$*_gc.log 2>logs/$*_gc.log

##################
# grab GTF files #
##################
gtfs/%.gtf:
	$(eval LOCAL_GTF_NAME=$(shell awk -F '\t' '$$1 == "$*" {print $$7}' names_experiments_table.txt))
	
	wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/gtf/$*/$(LOCAL_GTF_NAME) -O gtfs/$*.gtf.gz
	gunzip gtfs/$*.gtf.gz

###############
# find motifs #
###############
casx_sites_out/%_casx_sites.csv: fasta/%.fa
	motif_scraper --file_buffer --motif $(SEARCH_MOTIF) --cores $(MULTI_CPUS) --outputFile $@ $< 1>logs/$*_casx_search.log 2>&1
	
casx_sites_out/%_casx_sites.csv.gz: casx_sites_out/%_casx_sites.csv
	pigz --best -p $(MULTI_CPUS) $<

##############################################
# R data analysis                            #
# makes input for markdown figure generation #
##############################################
$(ROUTPUT): $(CASX_SITES_OUT) $(GTFFILES)
	R --vanilla < mdplyr_casx_analysis.R 1>logs/R_initial.log 2>&1 && \
	pigz --best -p $(MULTI_CPUS) R_analysis/*_cutsPerGene.csv && \
	pigz --best -p $(MULTI_CPUS) R_analysis/*_annotated_Casx_sites.csv && \
	pigz --best -p $(MULTI_CPUS) R_analysis/*_targets.csv
	
###########################
# pull mouse target sites #
###########################
R_analysis/mouse_target_upset.csv: $(ROUTPUT)
	R --vanilla < pull_mouse_casx_targets.R 1>logs/R_extract_mouse_targets.log 2>&1

$(MOUSE_TARGETS): R_analysis/mouse_target_upset.csv
	gzip --best $<
	
####################
# hamming distance #
####################
$(HAMMING): $(MOUSE_TARGETS) mouse_hamming_distance.py
	python mouse_hamming_distance.py

#########################################
# download FASTA files for each species #
#########################################
fasta/%.fa:
	wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/$*/dna/*dna_sm.primary_assembly.fa.gz -O fasta/$*.fa.gz
	
	if [ ! -s fasta/$*.fa.gz ]; \
	then wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/$*/dna/*dna_sm.toplevel.fa.gz -O fasta/$*.fa.gz; \
	fi;
	gunzip fasta/$*.fa.gz
