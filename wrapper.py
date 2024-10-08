#!/usr/bin/env python3

"""
Wrapper for splicing_analysis.py
"""

### ---------------------------------------- ###

def parse_argv():
    
    """
    This function parses command line arguments
    """

    print(f'args = {argv}')
    
    # Genes of interest, comma separated (N.B. must be ensembl ids)
    genes_list = argv[argv.index('--genes') + 1].split(',')
    
    # GTF file path
    gtf_path = argv[argv.index('--gtf') + 1]

    # Bam files, comma separated
    bam_paths = argv[argv.index('--bam') + 1].split(',')

    # Sample relationships (one must be "proband")
    # Must be in the same order as the bam files
    relationships = argv[argv.index('--relationships') + 1].split(',')

    # Filter junctions with poor coverage?
    if '--low_coverage_filter' in argv:
        
        kneedle_filter = True
    
    else:
        
        kneedle_filter = False
    
    # Trio VCF file (optional)
    # If provided, it is assumed that samples columns are in the same order as "relations"
    if '--vcf' in argv:
    
        vcf_path = argv[argv.index('--vcf') + 1]
    
    else:
        
        vcf_path = ''

    # Tab-separated file with the location of variants of interest
    if '--interesting_variants' in argv:
    
        interesting_variants = argv[argv.index('--interesting_variants') + 1]
    
    else:
        
        interesting_variants = ''
    
    return genes_list, gtf_path, bam_paths, relationships, kneedle_filter, vcf_path, interesting_variants

### ------------------MAIN------------------ ###

import splicing_analysis as sa

from sys import argv

### Parse args

genes_list, gtf_path, bam_paths, relationships, kneedle_filter, vcf_path, interesting_variants = parse_argv()

### Run analysis

analysis = sa.splice_check(genes_list, gtf_path, bam_paths, relationships, vcf_path, interesting_variants)

analysis.parse_gene_reads()

analysis.interpolate_junctions(kneedle_filter)
