#!/usr/bin/env python3

"""
"""

### ---------------------------------------- ###

class splice_check:

    """
    Class for ...

    Attributes
    ----------

    Methods
    -------
    """

    def __init__(self, genes, gtf_path, bam_paths, samples, trio_vcf_path='', variants_of_significance=''):

        # List of genes of interest
        self.genes = genes

        # Exon coordinates of genes of interest
        self.exon_coords = self.parse_gtf(gtf_path, genes)

        # Sample IDs
        self.sample_ids = samples

        # Paths of bam files
        self.bam_files = {s : bp for s,bp in zip(samples, bam_paths)}

        # Variants useful for distinguishing maternal and paternal reads (if a trio vcf is provided)
        if trio_vcf_path:

            genes_coords = {gene : (self.exon_coords.loc[self.exon_coords.gene_id == gene, 'contig'].values[0],
                            self.exon_coords.loc[self.exon_coords.gene_id == gene, 'start'].min(),
                            self.exon_coords.loc[self.exon_coords.gene_id == gene, 'stop'].max())
                           for gene in set(self.exon_coords.gene_id)}
    
            self.genes_coords = pd.DataFrame(genes_coords, index=['contig', 'start', 'stop']).T
            
            self.genes_vars = self.parse_vcf(trio_vcf_path)

        else:
            
            genes_coords = {gene : (self.exon_coords.loc[self.exon_coords.gene_id == gene, 'contig'].values[0],
                                    self.exon_coords.loc[self.exon_coords.gene_id == gene, 'start'].min(),
                                    self.exon_coords.loc[self.exon_coords.gene_id == gene, 'stop'].max())
                            for gene in set(self.exon_coords.gene_id)}
            
            self.genes_coords = pd.DataFrame(genes_coords, index=['contig', 'start', 'stop']).T
            
            self.genes_vars = {gene : pd.DataFrame([], columns=['contig', 'start', 'stop', 'allele', 'origin', 'var_type']) for gene in set(self.exon_coords.gene_id)}

        # Variants of interest to be visualized for each gene (if any)
        self.variants_of_significance = pd.read_csv(variants_of_significance, sep='\t')

    ### ------------------------------------ ###
    ### INPUT PARSING                        ###
    ### ------------------------------------ ###

    def parse_vcf(self, vcf_path):
        
        """
        A trio vcf file is parsed, keeping only variants of interest that can help distinguish between maternal and paternal
        """
        
        ### First round of screening: only keep variants in contigs of interest
        
        print('### Screening variants, first pass...')
        first_pass_vars = []
        tot_vars, useful_vars = 0, 0
        with gzip.open(vcf_path) as vcf:
            
            for line in vcf:
                
                line = line.decode('utf8')
                
                if line.startswith('##'):
                    
                    continue
                
                elif line.startswith('#CHROM'):
                    
                    colnames = line.replace('\n', '').split('\t')
                    
                    proband_id, parent_1_id, parent_2_id = colnames[-3:]
                
                else:
                    
                    tot_vars += 1
                    # Verbose counting
                    if tot_vars % 100000 == 0 and tot_vars > 0:
                        
                        print(f'Parsed {tot_vars} variants')
                    
                    contig, pos_start, _, ref, alt, _, _, _, sample_format, proband, parent_1, parent_2 = line.replace('\n', '').split('\t')
                    pos_start = int(pos_start)
                    pos_end = pos_start + len(ref)
                    
                    # Only keep variants in 
                    if contig not in self.genes_coords.contig.values:
                        
                        continue
                    
                    else:
                        
                        useful_vars += 1
                        first_pass_vars.append([contig, pos_start, pos_end, ref, alt, sample_format, proband, parent_1, parent_2])
            
            first_pass_vars = pd.DataFrame(first_pass_vars)
            first_pass_vars.columns = ['contig', 'pos_start', 'pos_end', 'ref', 'alt', 'sample_format', 'proband', 'parent_1', 'parent_2']
            
            print(f'{useful_vars} / {tot_vars} ({100 * useful_vars / (tot_vars)} %) variants were kept from first pass')
                
            ### Second round of screening: only keep variants useful for the genes of interest
            
            print('### Screening variants, second pass...')
            gene_vars = {}
            for gene, (contig, start, stop) in self.genes_coords.iterrows():
                
                # Subset variants
                vars_subset = first_pass_vars.loc[(first_pass_vars.contig == contig) &
                                                  (first_pass_vars.pos_start >= start) &
                                                  (first_pass_vars.pos_end <= stop),]
                
                # Check if variants can be used to discriminate maternal and paternal genes
                gene_vars[gene] = []
                for _, vs in vars_subset.iterrows():
                    
                    gene_vars[gene].extend(self.parse_variant(vs))
                
                gene_vars[gene] = pd.DataFrame(gene_vars[gene])
                gene_vars[gene].origin = [self.sample_ids[1] if orgn == 'parent_1' else self.sample_ids[2] if orgn == 'parent_2' else '' for orgn in gene_vars[gene].origin]
            
                print(f'Found { gene_vars[gene].shape[0]} good variants for gene {gene}')
        
        return gene_vars

    ### ------------------------------------ ###

    @staticmethod
    def parse_gtf(path, genes=[]):
        
        """
        This function parses a GTF file and returns the exon coordinates of genes of interest
        """
        
        print('### Parsing GTF file')
        
        # Load gtf annotation
        gtf = pd.read_csv(path, sep='\t', comment='#', header=None)
        gtf.columns = ['contig', 'source', 'feature', 'start', 'stop', 'score', 'strand', 'frame', 'attributes']
        
        # Keep only exons
        gtf = gtf.loc[gtf.feature == 'exon']
        
        # Extract gene ids
        field_start = 'gene_id "'
        field_stop = '"'
        gene_ids = [attr[attr.index(field_start) + len(field_start) :
                         attr[attr.index(field_start) + len(field_start) :].index(field_stop) + attr.index(field_start) + len(field_start)]
                    for attr in gtf.attributes]
        
        gtf = gtf.assign(gene_id = gene_ids)
        
        # Extract exon ids
        #field_start = 'exon_number "'
        #field_stop = '"'
        #exon_ids = [attr[attr.index(field_start) + len(field_start) :
        #                 attr[attr.index(field_start) + len(field_start) :].index(field_stop) + attr.index(field_start) + len(field_start)]
        #            for attr in gtf.attributes]
        exon_ids = np.arange(1, gtf.shape[0] + 1, 1)
        
        gtf = gtf.assign(exon_id = exon_ids)
        
        # Simplify table
        cols_to_keep = ['gene_id', 'exon_id', 'contig', 'start', 'stop']
        exon_coords = gtf.loc[:, cols_to_keep].copy()
        
        # Keep only desired genes, if specified
        if len(genes):
            
            exon_coords = exon_coords.loc[(exon_coords.gene_id.isin(genes)) |
                                          (exon_coords.gene_id.isin(genes)),]
        
        # Drop duplicate exons
        exon_coords.drop_duplicates(subset=['gene_id', 'contig', 'start', 'stop'], keep='first', inplace=True)
        
        return exon_coords

    ### ------------------------------------ ###

    @staticmethod
    def parse_variant(vs):
        
        """
        Parse variants and keep only those that can help distinguish between paternal- and maternal-derived reads
        """

        contig, pos_start, pos_end, ref, alt, sample_format, proband, parent_1, parent_2 = vs.values
        
        alleles = [ref] + alt.split(',')
        max_allele_len = max([len(a) for a in alleles])
        alleles_type = ['ref'] + ['snp' if (len(a) == len(ref)) & (len(ref) == 1) else 'indel' if len(a) < len(ref) else 'insertion' for a in alleles[1:]]
        alleles = [a.replace('*', '-') + '-' * (max_allele_len - len(a)) for a in alleles]
        
        genotype_index = sample_format.split(':').index('GT')
           
        proband_genotype, parent_1_genotype, parent_2_genotype = [sample.split(':')[genotype_index] for sample in [proband, parent_1, parent_2]]
        
        if '.' in f'{proband_genotype}/{parent_1_genotype}/{parent_2_genotype}':
            
            return []
        
        else:
            
            # Find alleles that allow to discriminate between paternal and maternal
            good_alleles = []
            for allele in set(proband_genotype.split('/')):
            
                if allele in parent_1_genotype and allele not in parent_2_genotype: # Parent 1 derived
                    
                    origin = 'parent_1'
                
                elif allele not in parent_1_genotype and allele in parent_2_genotype: # Parent 2 derived
                    
                    origin = 'parent_2'
                
                else:
                    
                    continue
                
                allele_seq = alleles[int(allele)]
                var_type = alleles_type[int(allele)]
                
                var = pd.Series({'contig' : contig,
                                 'start' : pos_start,
                                 'stop' : pos_end,
                                 'allele' : allele_seq,
                                 'origin' : origin,
                                 'var_type' : var_type})
                
                good_alleles.append(var)
                
            return good_alleles

    ### ------------------------------------ ###
    ### PARSE READS                          ###
    ### ------------------------------------ ###

    def parse_gene_reads(self):

        """
        Wrapper for parse_sample_reads
        """

        # Init results dicts
        self.junctions_dict = {}
        self.coverages_dict = {}
        self.reads_stats = {}
        
        # Init colors for sashimi plots
        colors = ['red', 'green', 'deepskyblue']

        # Process genes
        for gene, variants in self.genes_vars.items():
            
            print(f'### Parsing reads for {gene}')

            self.junctions_dict[gene] = {}
            self.coverages_dict[gene] = {}
            self.reads_stats[gene] = {}

            for n,(sample,bam_path) in enumerate(self.bam_files.items()):
                
                print(f'Processing {sample}')

                # Create output directory
                output_dir = f'{gene}/{sample}'

                if not exists(output_dir):
        
                    makedirs(output_dir)

                # Parse reads
                if sample == 'proband':
                    
                    read_stats, exon_coords, junctions, coverage = self.parse_sample_reads(gene, variants, sample, bam_path, colors[n])
                
                else:
                    
                    read_stats, exon_coords, junctions, coverage = self.parse_sample_reads(gene, pd.DataFrame([], columns=['contig', 'start', 'stop', 'allele', 'origin', 'var_type']), sample, bam_path, colors[n])

                self.junctions_dict[gene][sample] = junctions['unknown'].copy()
                self.coverages_dict[gene][sample] = coverage['unknown'].copy()
                self.reads_stats[gene][sample] = read_stats.copy()

            # Save stats to file
            pd.DataFrame(self.reads_stats[gene]).to_csv(f'{gene}/reads_origin_stats.tsv', sep='\t', index=True)

            # Sashimi plot for all samples together
            self.plot_all_sashimi(exon_coords, self.coverages_dict[gene], self.junctions_dict[gene], clean_junctions=False, out_name=f'{gene}/sashimi_all_samples.png')
            self.plot_all_sashimi(exon_coords, self.coverages_dict[gene], self.junctions_dict[gene], clean_junctions=True, out_name=f'{gene}/sashimi-cleaned_all_samples.png')

    ### ------------------------------------ ###

    def parse_sample_reads(self, gene, variants, sample_name, bam_path, color='red'):

        """
        Parsing aligned reads to find gene coverage and splicing events
        """

        # Open bam file
        alignments = pysam.AlignmentFile(bam_path, 'rb')

        # Get total number of mapped reads
        tot_sample_reads = sum([int(field.replace('mapped=', '').replace(',', ''))
                                for idxs in alignments.get_index_statistics()
                                for field in str(idxs).split(' ')
                                if field.startswith('mapped=')])

        # Get exons for specific gene
        gene_exon_coords = self.exon_coords.loc[self.exon_coords.gene_id == gene,].copy()

        # Save to file
        gene_exon_coords.to_csv(f'{gene}/{sample_name}/exons.tsv', sep='\t', index=False)
        
        # Fetch reads
        gene_contig, gene_start, gene_end = self.genes_coords.loc[gene,].values
        reads_of_interest = alignments.fetch(contig=gene_contig, start=gene_start - 1, stop=gene_end)
            
        # Parse reads and get parental origin (if vcf was provided)
        reads_info = {}
        for read in reads_of_interest:
                
            read_name, read_origin, ref_to_read_pos = self.get_read_origin(read, variants.copy())
                
            reads_info[read_name] = [read_origin, ref_to_read_pos]

        # Get number of reads mapping to the gene, then compute a normalization factor for read counts
        tot_gene_reads = len(reads_info)
        norm_factor = 1 / (tot_gene_reads * 1e6 / tot_sample_reads)
        
        # Collect stats based on read origin (if vcf was provided)
        if variants.shape[0]:
            
            origin_categories = ['unknown', 'possible_crossover', self.sample_ids[1], self.sample_ids[2]]
            read_category_stats = {category : sum([1 for val in reads_info.values() if val[0] == category]) for category in origin_categories}

        else:
            
            origin_categories = ['unknown']
            read_category_stats = {'unknown' : len(reads_info)}

        # For each category, format data and do a sashimi plot
        sample_gene_junctions, sample_gene_coverage = {}, {}
        for category in origin_categories:
            
            if category == 'unknown':
                
                # Get all the reads
                reads_covs = [ri[1] for ri in reads_info.values()]
            
            else:
                
                reads_covs = [ri[1] for ri in reads_info.values() if ri[0] == category]
            
            if not len(reads_covs):
                
                continue
            
            # Get coverage and splicing events
            gene_cov, unique_junctions = self.get_coverage_and_splicing(reads_covs, gene_exon_coords, norm_factor=norm_factor, out_name=f'{gene}/{sample_name}/junctions_origin={category}.tsv')
            
            # Store junctions and coverage data
            sample_gene_junctions[category] = unique_junctions.copy()
            sample_gene_coverage[category] = gene_cov.copy()

            # Sashimi plot
            self.plot_sample_sashimi(gene_exon_coords.copy(), gene_cov.copy(), unique_junctions.copy(), interval=[], out_name=f'{gene}/{sample_name}/sashimi_origin={category}.png', color=color)
            
            # Remove junctions not mapping to exons
            unique_junctions = unique_junctions.loc[(unique_junctions.left_exon != '') & (unique_junctions.right_exon != ''),]
            
            # Sashimi plot
            self.plot_sample_sashimi(gene_exon_coords.copy(), gene_cov.copy(), unique_junctions.copy(), interval=[], out_name=f'{gene}/{sample_name}/sashimi-cleaned_origin={category}.png', color=color)

        return read_category_stats, gene_exon_coords, sample_gene_junctions, sample_gene_coverage

    ### ------------------------------------ ###

    @staticmethod
    def get_coverage_and_splicing(r_covs, ex_coords, norm_factor=1, out_name='junctions.tsv'):

        """
        Compute gene coverage and detect splicing events
        """
    
        # Compute coverage
        gene_cov = np.concatenate(r_covs).astype(float)
        gene_cov = pd.DataFrame(gene_cov[(gene_cov[:, 0] != -1) & (gene_cov[:, 1] != -1), 1], columns=['pos'])
        gene_cov = gene_cov.groupby(by='pos').size()
        gene_cov = np.array([gene_cov.index, gene_cov.values])
        
        # Trim to ex_coords interval
        gene_start, gene_stop = ex_coords.start.min(), ex_coords.stop.max()
        gene_cov = gene_cov[:, (gene_cov[0,] >= gene_start) & (gene_cov[0,] <= gene_stop)]
        
        # Find splicing events as jumps in reads mapping
        junctions = []
        junctions_reads_proof = []
        for n,rc in enumerate(r_covs):
            
            # Remove -1 values in reads_positions
            rc_sub = rc[(rc[:,0] != -1) & (rc[:,1] != -1),]
            
            # Find junctions as regions where the reference position jumps
            diff = np.diff(rc_sub[:, 1])
            j_idx = np.where(diff > 1)[0]
            
            if len(j_idx):
                
                for i in j_idx:
                    
                    start, stop = rc_sub[i, 1], rc_sub[i + 1, 1]
                    
                    junctions.append([start, stop])
                    junctions_reads_proof.append(n)
            
        # Find unique junctions and the number of reads that support it
        unique_junctions = pd.DataFrame(junctions, columns=['start', 'stop']).groupby(by=['start', 'stop']).size().to_frame('reads_count').reset_index()
        
        # Normalize counts
        if norm_factor != 1:
            
            unique_junctions['norm_reads_count'] = unique_junctions['reads_count'].values * norm_factor
        
        # Annotate with possible exons
        exons_left, exons_right = [], []
        for _,(start, stop, *_) in unique_junctions.iterrows():
            
            left = ','.join(ex_coords.loc[ex_coords.stop == start, 'exon_id'].values.astype(str))
            right = ','.join(ex_coords.loc[ex_coords.start == stop, 'exon_id'].values.astype(str))
            
            exons_left.append(left)
            exons_right.append(right)
        
        unique_junctions['left_exon'] = exons_left
        unique_junctions['right_exon'] = exons_right
        
        # Save to file
        unique_junctions.to_csv(out_name, sep='\t', index=False, header=True)
        
        return gene_cov, unique_junctions

    ### ------------------------------------ ###

    @staticmethod
    def get_read_origin(read, var_data):

        """
        Give variants that can distinguish between paternal and maternal, reads are classified based on their origin
        """
        
        # Get read id
        read_name = read.query_name
        
        # Map reference and read positions
        ref_to_read_pos = np.array(read.get_aligned_pairs())
        ref_to_read_pos[ref_to_read_pos[:, 0] == None, 0] = -1
        ref_to_read_pos[ref_to_read_pos[:, 1] == None, 1] = -2
        ref_to_read_pos[:, 1] += 1
        
        # Prune var_data
        var_data = var_data.loc[(var_data.start >= ref_to_read_pos[:, 1].min()) & (var_data.stop <= ref_to_read_pos[:, 1].max()),]
        
        # Determine origin based on variants the read encompasses
        read_origin = []
        for _,vd in var_data.iterrows():
        
            if vd.start not in ref_to_read_pos[:, 1] or vd.stop - 1 not in ref_to_read_pos[:, 1]:
                
                continue
            
            else:
                
                target_allele = vd.allele.replace('-', '')
                
                # Find positions on the read that map to the coordinates of the variant
                read_pos_start = np.where(ref_to_read_pos[:, 1] == vd.start)[0][0]
                read_pos_stop = np.where(ref_to_read_pos[:, 1] == vd.stop - 1)[0][0]
                read_pos = ref_to_read_pos[read_pos_start : read_pos_stop + 1, 0]
                
                # Assemble variant from sequence
                read_var = ''.join([read.query_sequence[rpos] if rpos != -1 else '' for rpos in read_pos])
                
                if read_var == target_allele:
                    
                    origin = vd.origin
                    
                    # If variant is '*', then the positions leading and following the variant should not be empty, or it may simply be a splicing site
                    if not len(target_allele) and -1 in ref_to_read_pos[[read_pos_start - 1, read_pos_stop], 0]:
                        
                        origin = ''
                
                else:
                    
                    origin = ''
                
                if origin != '':
                    
                    read_origin.append(origin)
        
        if not len(read_origin):
            
            read_origin = 'unknown'
        
        elif len(set(read_origin)) == 1:
            
            read_origin = read_origin[0]
            
        elif len(set(read_origin)) == 2:
            
            read_origin = 'possible_crossover'
        
        else:
            
            read_origin = 'unknown'
        
        return read_name, read_origin, ref_to_read_pos

    ### ------------------------------------ ###
    ### INTERPOLATE JUNCTIONS DATA           ###
    ### ------------------------------------ ###
    
    def interpolate_junctions(self, kneedle_filter=True):
        
        """
        Wrapper for interpolate_gene_junctions
        """
        
        # Init results dicts
        self.family_junctions_dict = {}
        
        # Process genes
        for gene in self.genes:
            
            # Interpolate data
            print(f'### Interpolating junctions for {gene}')
            
            gene_family_junctions = self.interpolate_gene_junctions(gene, kneedle_filter)
            
            self.family_junctions_dict[gene] = gene_family_junctions.copy()
            
            # Get exons for specific gene
            gene_exon_coords = self.exon_coords.loc[self.exon_coords.gene_id == gene,].copy()
            
            # Only keep junctions with flags for plotting
            gene_family_junctions = gene_family_junctions.loc[(gene_family_junctions == 'FLAG').sum(axis=1) != 0,]
            
            # Plot data
            self.plot_sashimi_flags(gene_exon_coords.copy(), gene_family_junctions.copy(), self.variants_of_significance.copy(), out_name=f'{gene}/proband_flagged_junctions.png')
    
    ### ------------------------------------ ###
    
    def interpolate_gene_junctions(self, gene, kneedle_filter=True):

        """
        Interpolating gene splicing junctions detected in the trio
        """
        
        # For debug purposes
        try:
    
            del all_data
            
        except:
            
            pass
        
        # Merge data
        for sample,data in self.junctions_dict[gene].items():
            
            try:
                
                all_data = pd.merge(all_data, data,
                                    how='outer', on=['start', 'stop', 'left_exon', 'right_exon'])
                
                all_data.columns = all_data.columns.to_list()[:-2] + [f'reads_count_{sample}', f'norm_reads_count_{sample}']
            
            except:
                
                all_data = data.loc[:,['start', 'stop', 'left_exon', 'right_exon', 'reads_count', 'norm_reads_count']].copy()
                all_data.columns = ['start', 'stop', 'left_exon', 'right_exon', f'reads_count_{sample}', f'norm_reads_count_{sample}']
        
        # Fill missing count values
        all_data[all_data.columns[4:]] = all_data[all_data.columns[4:]].fillna(value=0)
        
        ### Filter junctions
        
        if kneedle_filter:
        
            thresholds = {s : self.kneedle(all_data[f'norm_reads_count_{s}'].values)[1] for s in self.sample_ids}
            row_filter = (all_data[[f'norm_reads_count_{s}'
                                    for s in self.sample_ids]].values < list(thresholds.values())).sum(axis=1) < len(self.sample_ids)
        
            all_data = all_data.loc[row_filter,]
        
        else:
            
            thresholds = {s : 0 for s in self.sample_ids}
        
        ### Reset index
        
        all_data = all_data.reset_index(drop=True)
        
        ### Tag junction
        
        # Find junctions unique to proband
        new_col_name = 'proband_unique'
        new_col_data = np.array(['' for _ in range(all_data.shape[0])], dtype='<U64')
        row_filter = ((all_data['norm_reads_count_proband'] > thresholds['proband']) &
                      ((all_data[[f'norm_reads_count_{s}' for s in self.sample_ids if s != 'proband']] <= list(thresholds.values())[1:]).sum(axis=1) == (len(self.sample_ids) - 1)))
        new_col_data[row_filter] = 'FLAG'
        all_data[new_col_name] = new_col_data
        
        # Find junctions missing in proband
        new_col_name = 'proband_missing'
        new_col_data = np.array(['' for _ in range(all_data.shape[0])], dtype='<U64')
        row_filter = ((all_data['norm_reads_count_proband'] <= thresholds['proband']) &
                      ((all_data[[f'norm_reads_count_{s}' for s in self.sample_ids if s != 'proband']] > list(thresholds.values())[1:]).sum(axis=1) == (len(self.sample_ids) - 1)))
        new_col_data[row_filter] = 'FLAG'
        all_data[new_col_name] = new_col_data
        
        # Find junctions that break correlation between samples
        for n,s1 in enumerate(self.sample_ids):
            
            for s2 in self.sample_ids[n + 1:]:
                
                new_col_name = f'{s1}_vs_{s2}'
                new_col_data = np.array(['' for _ in range(all_data.shape[0])], dtype='<U64')
        
                junctions_subset = all_data.loc[(all_data[f'norm_reads_count_{s1}'] > 0) & (all_data[f'norm_reads_count_{s2}'] > 0),]
                
                x_val = junctions_subset[f'norm_reads_count_{s1}'].values
                y_val = junctions_subset[f'norm_reads_count_{s2}'].values
                
                outliers_idx, r, p = self.find_outliers(x_val, y_val, s1, s2, f'{gene}/{s1}_{s2}_linregress.png')
                outliers_idx = junctions_subset.index[outliers_idx]
                
                new_col_data[outliers_idx] = 'FLAG'
                all_data[new_col_name] = new_col_data
        
        ### Write to file
        
        all_data.to_csv(f'{gene}/family_junctions.tsv', sep='\t', index=False)
        
        return all_data
    
    ### ------------------------------------ ###
    
    @staticmethod
    def find_outliers(x_val, y_val, x_name, y_name, diagnostic_plot_name='linregress.png'):

        """
        Assuming that the proband and a relative use mostly the same splicing junctions in a
        similar way, then the normalized counts for shared junctions should correlate
        A linear regression is fit and outliers are defined as being > 2 standard deviations from
        the fit
        """
        
        # Stats
        slope, intercept, r, p, std_err = linregress(x_val, y_val, alternative='two-sided')
        std = std_err * (len(x_val))**0.5
        y_predict = np.array([slope * x + intercept for x in x_val])
        
        # Define outliers based on > 2 std from best-fit line
        outliers_idx = np.where(abs(y_val - y_predict) > 2 * std)[0]
        
        if len(diagnostic_plot_name):
        
            # Plot data
            plt.figure(figsize=(5, 5))
            plt.title(f'r = {r:.3f}\np = {p:.2e}', loc='left')
            plt.xlabel(f'{x_name} counts')
            plt.ylabel(f'{y_name} counts')
            plt.scatter(x_val, y_val, color='lightgray', lw=0.5, edgecolor='black')
            plt.scatter(x_val[outliers_idx], y_val[outliers_idx], color='red', lw=0.5, edgecolor='black')
            plt.plot([x_val.min(), x_val.max()], [y_predict.min(), y_predict.max()], 'green')
            plt.plot([x_val.min(), x_val.max()], [y_predict.min() - 2 * std, y_predict.max() - 2 * std], 'black', linestyle='--')
            plt.plot([x_val.min(), x_val.max()], [y_predict.min() + 2 * std, y_predict.max() + 2 * std], 'black', linestyle='--')
            plt.savefig(diagnostic_plot_name, dpi=300)
            plt.close()
        
        return outliers_idx, r, p

    ### ------------------------------------ ###
    
    @staticmethod
    def kneedle(vector):
    
        """
        Kneedle to find threshold cutoff.
        """
        
        # Sort data
        vector = np.sort(vector)[::-1]
        
        # Find gradient and intercept
        x0, x1 = 0, len(vector)
        y0, y1 = max(vector), min(vector)
        gradient = (y1 - y0) / (x1 - x0)
        intercept = y0
        
        # Compute difference vector
        difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(vector)]
        
        # Find max of difference_vector and define cutoff
        cutoff_index = difference_vector.index(max(difference_vector))
        cutoff_value = vector[cutoff_index]
        
        return cutoff_index, cutoff_value
    
    ### ------------------------------------ ###
    ### PLOTTING                             ###
    ### ------------------------------------ ###
    
    @staticmethod
    def plot_all_sashimi(ex_coords, g_cov, u_junc, clean_junctions=False, out_name='sashimi.png'):
        
        """
        Sashimi and coverage plot for a specific gene and all samples
        """

        # Aestetic parameters
        intron_thickness = 1
        exon_thickness = 8
        splicing_thickness = 1
        min_height = 1
        max_height = 10
        colors = ['deepskyblue', 'green', 'red']
        
        gene_start, gene_stop = ex_coords.start.min(), ex_coords.stop.max()

        fig_width = max(10, abs(gene_stop - gene_start) / 10000)
        fig_height = 12
        
        # Clean junctions by removing those not mapping to exons
        if clean_junctions:
            
            u_junc = {sample : j.loc[(j.left_exon != '') & (j.right_exon != ''),].copy() for sample,j in u_junc.items()}
        
        # Add splicing size to u_junc, normalized
        juncs_plot_data = pd.DataFrame([[sample, start, stop, abs(stop - start)]
                                        for sample,j in u_junc.items()
                                        for _,(start, stop, count, *_) in j.iterrows()],
                                       columns=['sample', 'start', 'stop', 'splicing_size'])
        juncs_plot_data = juncs_plot_data.assign(splicing_size_norm = min_height + np.log2(juncs_plot_data.splicing_size / juncs_plot_data.splicing_size.min()))
        
        # Normalize coverage to 0-1 range, then set g_cov_norm max height to 2 * juncs_plot_data.splicing_size_norm.max()
        g_cov_max = max([c[1,].max() for c in g_cov.values()])
        g_cov_norm = {}
        for sample,cov in g_cov.items():
            
            g_cov_norm[sample] = cov.copy()
            
            g_cov_norm[sample][1,] = (g_cov_norm[sample][1,] / g_cov_max)
            g_cov_norm[sample][1,] = g_cov_norm[sample][1,] * (2 * juncs_plot_data.splicing_size_norm.max()) if len(juncs_plot_data) else max_height
        
        # Init plot
        plt.figure(figsize=(fig_width, fig_height))
        
        ### Gene structure
        
        # Plot gene length
        gene_structure_pos = -3
        plt.hlines(gene_structure_pos, gene_start, gene_stop, linestyles='solid', colors='black', linewidth=intron_thickness)
        
        # Plot exons
        for _, exon in ex_coords.iterrows():
            
            plt.hlines(gene_structure_pos, exon.start, exon.stop, linestyles='solid', colors='black', linewidth=exon_thickness)
        
        ### Coverage and splicing
        
        offset = 0
        
        ordered_samples = list(g_cov_norm.keys())[::-1]
        
        for n,(s) in enumerate(ordered_samples):
            
            c = g_cov_norm[s].copy()
            
            j = juncs_plot_data.loc[juncs_plot_data['sample'] == s,]
            
            offset += 3
            
            offset += j.loc[((j.start >= gene_start) & (j.start <= gene_stop)) |
                            ((j.stop >= gene_start) & (j.stop <= gene_stop)), 'splicing_size_norm'].max()
        
            plt.hlines(offset, gene_start, gene_stop, linestyles='solid', colors='black', linewidth=intron_thickness)
            
            # Plot coverage
            where_pos = [True] + list(np.diff(c[0,]) == 1) # Where to fill
            plt.fill_between(x=c[0,], y1=np.zeros(c.shape[1]) + offset, y2=c[1,] + offset, where=where_pos, alpha=1, color=colors[n], linewidth=0)
            
            # Plot splicing events
            for _,splicing in j.iterrows():
                
                intercept = - splicing.splicing_size_norm
                
                if splicing.start < splicing.stop:
                    
                    spline = CubicSpline(x=[splicing.start, (splicing.stop + splicing.start) / 2, splicing.stop],
                                         y=[0, intercept, 0])
                
                else:
                    
                    spline = CubicSpline(x=[splicing.stop, (splicing.stop + splicing.start) / 2, splicing.start],
                                         y=[0, intercept, 0])
                
                x = np.linspace(splicing.start, splicing.stop, 21)
                y = spline(x) + offset
                
                plt.plot(x, y, ls='-', lw=splicing_thickness, color='black', alpha=1)
            
            offset += c[1,].max()
        
        plt.xlabel(f'{ex_coords.contig.values[0]} (bp)')
        plt.xlim(gene_start, gene_stop)
        plt.yticks([])
        plt.title(None)
        plt.tight_layout()
        plt.box(False)
        
        plt.savefig(out_name, dpi=300)
        plt.close()
    
    ### ------------------------------------ ###

    @staticmethod
    def plot_sample_sashimi(ex_coords, g_cov, u_junc, interval=[], out_name='sashimi.png', color='red'):
        
        """
        Sashimi and coverage plot for a specific gene and individual sample
        """
    
        # Aestetic parameters
        intron_thickness = 1
        exon_thickness = 5
        splicing_thickness = 1
        min_height = 1
        max_height = 10
        
        if not len(interval):
            
            gene_start, gene_stop = ex_coords.start.min(), ex_coords.stop.max()
            
        else:
            
            gene_start, gene_stop = interval
        
        fig_width = max(10, abs(gene_stop - gene_start) / 10000)
        fig_height = 4
        
        # Add splicing size to u_junc, normalized
        juncs_plot_data = pd.DataFrame([[start, stop, abs(stop - start)]
                                        for _,(start, stop, count, *_) in u_junc.iterrows()],
                                       columns=['start', 'stop', 'splicing_size'])
        juncs_plot_data = juncs_plot_data.assign(splicing_size_norm = min_height + np.log2(juncs_plot_data.splicing_size / juncs_plot_data.splicing_size.min()))
        
        # Normalize coverage to 0-1 range
        g_cov[1,] = (g_cov[1,] / g_cov[1,].max())
        
        # Set g_cov max height to 2 * juncs_plot_data.splicing_size_norm.max()
        g_cov[1,] = g_cov[1,] * (2 * juncs_plot_data.splicing_size_norm.max()) if len(juncs_plot_data) else max_height
        
        # Init plot
        plt.figure(figsize=(fig_width, fig_height))
        
        ### Gene structure
        
        if len(juncs_plot_data):
            
            gene_structure_pos = - juncs_plot_data.loc[((juncs_plot_data.start >= gene_start) & (juncs_plot_data.start <= gene_stop)) |
                                                       ((juncs_plot_data.stop >= gene_start) & (juncs_plot_data.stop <= gene_stop)), 'splicing_size_norm'].max() - 2
        
        else:
            
            gene_structure_pos = - 1
        
        # Plot gene length
        plt.hlines(gene_structure_pos, gene_start, gene_stop, linestyles='solid', colors='black', linewidth=intron_thickness)
        
        # Plot exons
        for _, exon in ex_coords.iterrows():
            
            plt.hlines(gene_structure_pos, exon.start, exon.stop, linestyles='solid', colors='black', linewidth=exon_thickness)
        
        ### Coverage and splicing
        
        plt.hlines(0, gene_start, gene_stop, linestyles='solid', colors='black', linewidth=intron_thickness)
        
        # Plot coverage
        where_pos = [True] + list(np.diff(g_cov[0,]) == 1) # Where to fill
        plt.fill_between(x=g_cov[0,], y1=np.zeros(g_cov.shape[1]), y2=g_cov[1,], where=where_pos, alpha=1, color=color, linewidth=0)
        
        # Plot splicing events
        for _,splicing in juncs_plot_data.iterrows():
            
            intercept = - splicing.splicing_size_norm
            
            if splicing.start < splicing.stop:
                
                spline = CubicSpline(x=[splicing.start, (splicing.stop + splicing.start) / 2, splicing.stop],
                                     y=[0, intercept, 0])
            
            else:
                
                spline = CubicSpline(x=[splicing.stop, (splicing.stop + splicing.start) / 2, splicing.start],
                                     y=[0, intercept, 0])
            
            x = np.linspace(splicing.start, splicing.stop, 21)
            y = spline(x)
            
            plt.plot(x, y, ls='-', lw=splicing_thickness, color='black', alpha=1)
        
        plt.xlabel(f'{ex_coords.contig.values[0]} (bp)')
        plt.xlim(gene_start, gene_stop)
        plt.yticks([])
        plt.title(None)
        plt.tight_layout()
        plt.box(False)
        
        plt.savefig(out_name, dpi=300)
        plt.close()

    ### ------------------------------------ ###
    
    @staticmethod
    def plot_sashimi_flags(ex_coords, juncs, interesting_vars, out_name='proband_flagged_junctions.png'):
        
        """
        Sashimi plot for the flagged junctions of a specific gene
        """
        
        # Aestetic parameters
        intron_thickness = 1
        exon_thickness = 5
        splicing_thickness = 1
        min_height = 1
        
        contig, gene_start, gene_stop = ex_coords.contig.values[0], ex_coords.start.min(), ex_coords.stop.max()
        
        interesting_vars = interesting_vars.loc[(interesting_vars.contig == contig) &
                                                (interesting_vars.pos >= gene_start) &
                                                (interesting_vars.pos <= gene_stop),]
        
        juncs = juncs.loc[((juncs.start >= gene_start) & (juncs.start <= gene_stop)) |
                          ((juncs.stop >= gene_start) & (juncs.stop <= gene_stop)),]
        
        fig_width = max(10, abs(gene_stop - gene_start) / 10000)
        fig_height = 4
        
        # Add splicing size to juncs, normalized
        juncs_plot_data = pd.DataFrame([[start, stop, abs(stop - start)] for _,(start, stop, count, *_) in juncs.iterrows()], columns=['start', 'stop', 'splicing_size'])
        juncs_plot_data = juncs_plot_data.assign(splicing_size_norm = min_height + np.log2(juncs_plot_data.splicing_size / juncs_plot_data.splicing_size.min()))
        
        # Init plot
        plt.figure(figsize=(fig_width, fig_height))
        
        ### Gene structure
        
        gene_structure_pos = - 1
        
        # Plot gene length
        plt.hlines(gene_structure_pos, gene_start, gene_stop, linestyles='solid', colors='black', linewidth=intron_thickness)
        
        # Plot exons
        for _, exon in ex_coords.iterrows():
            
            plt.hlines(gene_structure_pos, exon.start, exon.stop, linestyles='solid', colors='black', linewidth=exon_thickness)
        
        ### Coverage and splicing
        
        plt.hlines(0, gene_start, gene_stop, linestyles='solid', colors='black', linewidth=intron_thickness)
        
        # Plot splicing events
        for _,splicing in juncs_plot_data.iterrows():
            
            intercept = splicing.splicing_size_norm
            
            if splicing.start < splicing.stop:
                
                spline = CubicSpline(x=[splicing.start, (splicing.stop + splicing.start) / 2, splicing.stop],
                                     y=[0, intercept, 0])
            
            else:
                
                spline = CubicSpline(x=[splicing.stop, (splicing.stop + splicing.start) / 2, splicing.start],
                                     y=[0, intercept, 0])
            
            x = np.linspace(splicing.start, splicing.stop, 21)
            y = spline(x)
            
            plt.plot(x, y, ls='-', lw=splicing_thickness, color='black', alpha=1)
        
        ### FLAGS
        
        # Find columns with flags
        flags_cols = [col for col in juncs.columns[np.sort(np.unique(np.where(juncs == 'FLAG')[1]))] if 'proband' in col]
        
        # Find x and y axis positions as well as flag type
        base_y = gene_structure_pos = juncs_plot_data['splicing_size_norm'].max() + 2
        xvals, yvals, flags = [], [], []
        for _,j in juncs.iterrows():
            
            x = (j.start + j.stop) / 2
            for n,col in enumerate(flags_cols):
                
                if j[col] == 'FLAG':
                    
                    y = base_y + n
                    
                    xvals.append(x)
                    yvals.append(y)
                    flags.append(col)
        
        # Add candidate variant possibly causing splicing abnormalities
        for _,(_,pos,flag) in interesting_vars.iterrows():
            
            xvals.append(pos)
            yvals.append(base_y + len(flags_cols))
            flags.append(flag)
        
        # Convert to data frame
        flags_plot_data = pd.DataFrame(np.stack([xvals, yvals, flags], axis=-1), columns=['x', 'y', 'flags'])
        flags_plot_data = flags_plot_data.astype({'x' : float, 'y' : float, 'flags' : str})
        
        # Plotting flags
        scatter = sns.scatterplot(flags_plot_data, x='x', y='y', hue='flags', s=25, edgecolor='black', linewidth=1)
        sns.move_legend(scatter, loc='upper left', bbox_to_anchor=(1, 1), title='Flags')
        
        plt.xlabel(f'{ex_coords.contig.values[0]} (bp)')
        plt.xlim(gene_start, gene_stop)
        plt.ylabel(None)
        plt.yticks([])
        plt.title(None)
        plt.tight_layout()
        plt.box(False)
        
        plt.savefig(out_name, dpi=300)
        plt.close()

### ------------------MAIN------------------ ###

import gzip
import numpy as np
import pandas as pd
import pysam
import seaborn as sns

from matplotlib import pyplot as plt
from os import makedirs
from os.path import exists
from scipy.interpolate import CubicSpline
from scipy.stats import linregress
