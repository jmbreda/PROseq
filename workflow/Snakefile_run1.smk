
configfile: "config/config_run1.yaml"

wildcard_constraints:
    sample="|".join(config['Samples']),
    strand="|".join(config['Strands']),
    bin_size="|".join(config['Bin_size'])

rule all:
    input:
        expand(os.path.join(config['outfold_star'],"{sample}","Aligned.sortedByCoord.out.bam"), sample=config['Samples']),
        #expand(os.path.join(config['outfold_binned'],"{sample}","NormCoverage_3p_{strand}_bin{bin_size}bp.bw"), sample=config['Samples'], strand=config['Strands'], bin_size=config['Bin_size'])
        #expand(os.path.join(config['outfold_binned'],"NormCoverage_3p_bin{bin_size}bp_{chr}.csv"),bin_size=config['Bin_size'],chr=config['Chromosomes']['GRCm39']),
        #expand(os.path.join(config['outfold_phase_amp'],"overall_phase_amp_{bin_size}bp.csv"), bin_size=config['Bin_size']),
        #expand(os.path.join(config['outfold_phase_amp'],"bin_phase_amp_{bin_size}bp_{strand}.bb"), bin_size=config['Bin_size'], strand=config['Strands']),
        #expand(os.path.join(config['outfold_phase_amp'],"gene_phase_amp_{bin_size}bp.bb"),bin_size=config['Bin_size'])
        #expand(os.path.join(config['outfold_fig'],"space_time_fourier","Bins_{bin_size}bp_sigma_{sigma}Mb","{chr}_amp_pos_k.pdf"),bin_size=config['Bin_size'],sigma=config['Space_Time_Fourier']['sigma_Mb'],chr=config['Chromosomes']['GRCm39'])
        #expand(os.path.join(config['outfold_fig'],"space_time_fourier","Bins_{bin_size}bp_sigma_{sigma}Mb","chr17_amp_pos_k.pdf"),bin_size=config['Bin_size'],sigma=config['Space_Time_Fourier']['sigma_Mb']),
        #expand(os.path.join(config['outfold_kalman'],"ll_{win_size}_{d_win}_{chr}_{bin_size}bp.csv"),win_size=config['Kalman']['win_size'],d_win=config['Kalman']['d_win'],chr=config['Chromosomes']['GRCm39'],bin_size=config['Bin_size']),
        #expand(os.path.join(config['outfold_kalman'],"mu_{win_size}_{d_win}_{chr}_{bin_size}bp.npy"),win_size=config['Kalman']['win_size'],d_win=config['Kalman']['d_win'],chr=config['Chromosomes']['GRCm39'],bin_size=config['Bin_size']),
        #expand(os.path.join(config['outfold_kalman'],"sigma_{win_size}_{d_win}_{chr}_{bin_size}bp.npy"),win_size=config['Kalman']['win_size'],d_win=config['Kalman']['d_win'],chr=config['Chromosomes']['GRCm39'],bin_size=config['Bin_size']),
        #expand(os.path.join(config['outfold_kalman'],"ll_{win_size}_{d_win}_{chr}_{bin_size}bp.csv"),win_size=config['Kalman']['win_size'],d_win=config['Kalman']['d_win'],chr=config['Chromosomes']['Test'],bin_size=config['Bin_size']),
        #expand(os.path.join(config['outfold_kalman'],"mu_{win_size}_{d_win}_{chr}_{bin_size}bp.csv"),win_size=config['Kalman']['win_size'],d_win=config['Kalman']['d_win'],chr=config['Chromosomes']['Test'],bin_size=config['Bin_size']),
        #expand(os.path.join(config['outfold_kalman'],"sigma_{win_size}_{d_win}_{chr}_{bin_size}bp.csv"),win_size=config['Kalman']['win_size'],d_win=config['Kalman']['d_win'],chr=config['Chromosomes']['Test'],bin_size=config['Bin_size'])
        #expand(os.path.join(config['outfold_fig'],"kalman","ll_{win_size}_{d_win}_{chr}_{bin_size}bp.pdf"),win_size=config['Kalman']['win_size'],d_win=config['Kalman']['d_win'],chr=config['Chromosomes']['GRCm39'],bin_size=config['Bin_size'])
        #expand( os.path.join(config['outfold_binned'],"Noise_model_parameters_{bin_size}bp.csv"), bin_size=config['Bin_size'] )
        #expand(os.path.join(config['outfold_kalman'],"Gene_Q_1e-3","Kalman_Smoothing_{bin_size}bp.hdf5"), bin_size=config['Bin_size']),
        expand(os.path.join(config['outfold_fig'],"kalman","Gene_Q_1e-3","scatter_Kmax_R2_LL_strand_{bin_size}bp.pdf"), bin_size=config['Bin_size'])

##-------------------------------------------------------------------##
##   Download gtf, cDNA (coding transcripts), and DNA from Gencode   ##
##-------------------------------------------------------------------##

rule download_ensembl_gtf:
    output: 
        os.path.join(config['outfold_genome'],config['Species'],"gene_annotation.gtf")
    params: 
        url = lambda wildcards: config['URLs']['GRCm39']['gtf']
    shell:
        """
        wget -O {output}.gz {params.url}
        gunzip {output}.gz
        """

rule download_ensembl_cdna:
    output: 
        os.path.join(config['outfold_genome'],config['Species'],"cdna.fa")
    params: 
        url = lambda wildcards: config['URLs']['GRCm39']['cdna']
    shell:
        """
        wget -O {output}.gz {params.url}
        gunzip {output}.gz
        """

rule download_ensembl_genome:
    output: 
        os.path.join(config['outfold_genome'],config['Species'],"genome.fa")
    params:
        url = lambda wildcards: config['URLs']['GRCm39']['genome']
    shell:
        """
        wget -O {output}.gz {params.url}
        gunzip {output}.gz
        """

##----------------------------------------##
##  Generate STAR genome index and align  ##
##----------------------------------------##

rule star_genome_index:
    input:
        genome = os.path.join(config['outfold_genome'],config['Species'],"genome.fa"),
        gtf = os.path.join(config['outfold_genome'],config['Species'],"gene_annotation.gtf")
    output:
        genome_dir = directory(os.path.join(config['outfold_genome'],config['Species'],"star_index"))
    params:
        overhang = config['read_length']-1
    threads: 24
    shell:
        """
        module load gcc/11.3.0 star
        mkdir -p {output.genome_dir}
        STAR --runMode genomeGenerate \
             --genomeDir {output.genome_dir} \
             --genomeFastaFiles {input.genome} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.overhang} \
             --limitGenomeGenerateRAM 81000000000 \
             --runThreadN {threads}
        """

rule cutadapt:
    input:
        fq = os.path.join(config['outfold_fastq'],"{sample}.fastq")
    output:
        fq_trimmed = os.path.join(config['outfold_fastq'],"{sample}_trimmed.fastq"),
        report = os.path.join(config['outfold_fastq'],"{sample}_trimmed_report.txt")
    params:
        adapter = config['adapter']
    threads: 12
    shell:
        """
        source activate cutadapt
        cutadapt -a {params.adapter} -o {output.fq_trimmed} -j 12 {input.fq} > {output.report}
        """

rule star_align:
    input:
        fq = os.path.join(config['outfold_fastq'],"{sample}_trimmed.fastq"),
        gtf = os.path.join(config['outfold_genome'],config['Species'],"gene_annotation.gtf"),
        genome_dir = os.path.join(config['outfold_genome'],config['Species'],"star_index")
    output:
        os.path.join(config['outfold_star'],"{sample}","Aligned.sortedByCoord.out.bam")
    params:
        prefix=os.path.join(config['outfold_star'],"{sample}/"),
        outsamtype = "BAM SortedByCoordinate",
        limitram = 81000000000
    threads: 12
    shell:
        """
        module load gcc/11.3.0 star
        STAR --runMode alignReads \
             --readFilesIn {input.fq} \
             --readFilesCommand gunzip -c \
             --genomeDir {input.genome_dir} \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype {params.outsamtype} \
             --limitBAMsortRAM {params.limitram} \
             --sjdbGTFfile {input.gtf} \
             --runThreadN {threads}
        """

##--------------------------------------##
##  BAM files indexing                  ##
##--------------------------------------##

rule sam_index:
    input:
        bam=os.path.join(config['outfold_star'],"{sample}/Aligned.sortedByCoord.out.bam")
    output:
        bam_index=os.path.join(config['outfold_star'],"{sample}/Aligned.sortedByCoord.out.bam.bai")
    shell:
        """
        ml gcc/11.3.0 samtools
        samtools index {input}
        """

##--------------------------------------##
##  Coverage                            ##
##--------------------------------------##

rule chrom_size:
    input:
        genome=os.path.join(config['outfold_genome'],config['Species'],"genome.fa")
    output:
        chrom_size=os.path.join(config['outfold_genome'],config['Species'],"mm39.chrom.sizes")
    shell:
        """
        ml gcc/11.3.0 samtools
        samtools faidx {input.genome}
        cut -f1,2 {input.genome}.fai | awk -v bs=1000 '{{print $1"\t"$2 + bs - $2 % bs}}' > {output.chrom_size}
        """

# make sure chromosome length are multiples of bin size
rule chrom_size_bin_size:
    input:
        genome=os.path.join(config['outfold_genome'],config['Species'],"genome.fa")
    output:
        chrom_size=os.path.join(config['outfold_genome'],config['Species'],"mm39.chrom.sizes.{bin_size}bp")
    shell:
        """
        ml gcc/11.3.0 samtools
        samtools faidx {input.genome}
        cut -f1,2 {input.genome}.fai | awk -v bs={wildcards.bin_size} '{{print $1"\t"$2 + bs - $2 % bs}}' > {output.chrom_size}
        """

# for PRO-seq data, we only need the 3'-end of the reads. However the protocol repots reverse complement of the reads. So I take the 5'-end and flip the strand.
rule coverage_bedgraph:
    input:
        bam=os.path.join(config['outfold_star'],"{sample}","Aligned.sortedByCoord.out.bam")
    output:
        bg=os.path.join(config['outfold_counts'],"{sample}","coverage_3p_{strand}.bedgraph")
    params:
        strand = lambda wildcards: '-' if wildcards.strand == "forward" else '+' if wildcards.strand == "reverse" else ''
    shell:
        """
        ml gcc/11.3.0 samtools bedtools2
        samtools view -b {input.bam} | bedtools genomecov -ibam stdin -bg -strand {params.strand} -5 | grep "^\<chr" > {output.bg}
        """

rule get_total_count:
    input:
        bg=os.path.join(config['outfold_counts'],"{sample}","coverage_3p_{strand}.bedgraph")
    output:
        counts=os.path.join(config['outfold_counts'],"{sample}","total_counts_{strand}.txt")
    shell:
        """
        count=$(awk '{{sum+=$4}} END {{print sum}}' {input.bg})
        echo $count > {output.counts}
        """

rule normalize_sort_bedgraph:
    input:
        bg=os.path.join(config['outfold_counts'],"{sample}","coverage_3p_{strand}.bedgraph"),
        total_count=os.path.join(config['outfold_counts'],"{sample}","total_counts_{strand}.txt"),
        total_counts=expand(os.path.join(config['outfold_counts'],"{sample}","total_counts_{strand}.txt"),sample=config['Samples'], strand=config['Strands'])
    output:
        norm_bg=os.path.join(config['outfold_norm'],"{sample}","NormCoverage_3p_{strand}.bedgraph")
    shell:
        """
        ./scripts/normalize_bedgraph.sh {input.bg} {input.total_count} {output.norm_bg} {input.total_counts}
        """

rule norm_coverage_bw:
    input:
        norm_bg=os.path.join(config['outfold_norm'],"{sample}","NormCoverage_3p_{strand}.bedgraph"),
        chrom_size=os.path.join(config['outfold_genome'],config['Species'],"mm39.chrom.sizes")
    output:
        bw=os.path.join(config['outfold_norm'],"{sample}","NormCoverage_3p_{strand}.bw")
    shell:
        """
        bedGraphToBigWig {input.norm_bg} {input.chrom_size} {output.bw}
        """

rule bin_coverage:
    input:
        bw=os.path.join(config['outfold_norm'],"{sample}","NormCoverage_3p_{strand}.bw")
    output:
        bin_bw=os.path.join(config['outfold_binned'],"{sample}","NormCoverage_3p_{strand}_bin{bin_size}bp.bw")
    shell:
        """
        python scripts/bin_bigwig.py  -i {input.bw} -o {output.bin_bw} -b {wildcards.bin_size}
        """

rule make_tables:
    input:
        expand(os.path.join(config['outfold_binned'],"{sample}","NormCoverage_3p_{strand}_bin{{bin_size}}bp.bw"), sample=config['Samples'], strand=config['Strands'])
    output:
        table=os.path.join(config['outfold_binned'],"NormCoverage_3p_bin{bin_size}bp_{chr}.csv")
    params:
        bw_folder=config['outfold_binned']
    shell:
        """
        python scripts/make_expression_table.py --bin_size {wildcards.bin_size} \
                                                --chr {wildcards.chr} \
                                                --bw_folder {params.bw_folder} \
                                                --output {output.table}
        """


##--------------------------------------##
##  Get bin phase and amplitude         ##
##--------------------------------------##

rule get_overall_phase_amp:
    input:
        expand(os.path.join(config['outfold_binned'],"{sample}","NormCoverage_3p_{strand}_bin{{bin_size}}bp.bw"), sample=config['Samples'], strand=config['Strands']) # only for dependency
    output:
        table=os.path.join(config['outfold_phase_amp'],"overall_phase_amp_{bin_size}bp.csv")
    params:
        bw_folder=config['outfold_binned']
    shell:
        """
        python scripts/get_overall_phase_amp.py --bin_size {wildcards.bin_size} \
                                                --bw_folder {params.bw_folder} \
                                                --out_table {output.table}
        """ 
    
rule get_bin_phase_amp:
    input:
        overall_phase_amp=os.path.join(config['outfold_phase_amp'],"overall_phase_amp_{bin_size}bp.csv")
    output:
        table=os.path.join(config['outfold_phase_amp'],"bin_phase_amp_{bin_size}bp.csv")
    params:
        bw_folder=config['outfold_binned']
    threads: 12
    shell:
        """
        python scripts/get_bin_phase_amp.py --bin_size {wildcards.bin_size} \
                                            --bw_folder {params.bw_folder} \
                                            --overall_phase_amp_table {input.overall_phase_amp} \
                                            --out_table {output.table}
        """

rule get_bin_phase_amp_bed:
    input:
        table=os.path.join(config['outfold_phase_amp'],"bin_phase_amp_{bin_size}bp.csv")
    output:
        bed=os.path.join(config['outfold_phase_amp'],"bin_phase_amp_{bin_size}bp.bed")
    shell:
        """
        python scripts/get_bin_phase_amp_bed_track.py --in_table {input.table} \
                                                      --out_bed {output.bed}
        """

rule separate_stands:
    input:
        bed=os.path.join(config['outfold_phase_amp'],"bin_phase_amp_{bin_size}bp.bed")
    output:
        bed_forward=os.path.join(config['outfold_phase_amp'],"bin_phase_amp_{bin_size}bp_forward.bed"),
        bed_reverse=os.path.join(config['outfold_phase_amp'],"bin_phase_amp_{bin_size}bp_reverse.bed")
    shell:
        """
        awk '{{if ($6 == "+") print $0}}' {input.bed} > {output.bed_forward}
        awk '{{if ($6 == "-") print $0}}' {input.bed} > {output.bed_reverse}
        """

rule bin_bed_to_bigBed:
    input:
        bed=os.path.join(config['outfold_phase_amp'],"bin_phase_amp_{bin_size}bp_{strand}.bed"),
        chrom_size=os.path.join(config['outfold_genome'],config['Species'],"mm39.chrom.sizes.{bin_size}bp")
    output:
        bb=os.path.join(config['outfold_phase_amp'],"bin_phase_amp_{bin_size}bp_{strand}.bb")
    shell:
        """
        bedSort {input.bed} {input.bed}
        bedToBigBed {input.bed} {input.chrom_size} {output.bb}
        """

##--------------------------------------##
##  Get gene phase and amplitude        ##
##--------------------------------------##

rule get_gene_gtf:
    input:
        gtf=os.path.join(config['outfold_genome'],config['Species'],"gene_annotation.gtf")
    output:
        gtf_gene=os.path.join(config['outfold_genome'],config['Species'],"gene_protein_coding.gtf")
    shell:
        """
        ./scripts/get_gene_gtf.sh {input.gtf} {output.gtf_gene}
        """

rule get_gene_phase_amp:
    input:
        gtf=os.path.join(config['outfold_genome'],config['Species'],"gene_protein_coding.gtf"),
        overall_phase_amp=os.path.join(config['outfold_phase_amp'],"overall_phase_amp_{bin_size}bp.csv")
    output:
        table=os.path.join(config['outfold_phase_amp'],"gene_phase_amp_{bin_size}bp.csv")
    params:
        bw_folder=config['outfold_binned']
    threads: 12
    shell:
        """
        python scripts/get_gene_phase_amp.py --gtf {input.gtf} \
                                             --bin_size {wildcards.bin_size} \
                                             --bw_folder {params.bw_folder} \
                                             --overall_phase_amp_table {input.overall_phase_amp} \
                                             --out_table {output.table}
        """

rule get_gene_phase_amp_bed:
    input:
        table=os.path.join(config['outfold_phase_amp'],"gene_phase_amp_{bin_size}bp.csv")
    output:
        bed=os.path.join(config['outfold_phase_amp'],"gene_phase_amp_{bin_size}bp.bed"),
        fig=os.path.join(config['outfold_fig'],"gene_phase_amp_{bin_size}bp.pdf")
    shell:
        """
        python scripts/get_gene_phase_amp_bed_track.py --in_table {input.table} \
                                                       --out_bed {output.bed} \
                                                       --out_fig {output.fig}
        """

rule gene_bed_to_bigBed:
    input:
        bed=os.path.join(config['outfold_phase_amp'],"gene_phase_amp_{bin_size}bp.bed"),
        chrom_size=os.path.join(config['outfold_genome'],config['Species'],"mm39.chrom.sizes.{bin_size}bp")
    output:
        bb=os.path.join(config['outfold_phase_amp'],"gene_phase_amp_{bin_size}bp.bb")
    shell:
        """
        bedSort {input.bed} {input.bed}
        bedToBigBed {input.bed} {input.chrom_size} {output.bb}
        """

##--------------------------------------##
##  Time-space Fourier                  ##
##--------------------------------------##


rule time_space_fourier_transform:
    input:
        overall_phase_amp=os.path.join(config['outfold_phase_amp'],"overall_phase_amp_{bin_size}bp.csv")
    output:
        amplitude=os.path.join(config['outfold_space_time_fourier'],"Bins_{bin_size}bp_sigma_{sigma}Mb","{chr}_amp_pos_k.csv"),
        phase=os.path.join(config['outfold_space_time_fourier'],"Bins_{bin_size}bp_sigma_{sigma}Mb","{chr}_phase_pos_k.csv"),
        mu=os.path.join(config['outfold_space_time_fourier'],"Bins_{bin_size}bp_sigma_{sigma}Mb","{chr}_mu_pos.csv")
    params:
        bw_folder=config['outfold_binned'],
    threads: 24
    shell:
        """
        python scripts/time_space_Fourrier_transform.py --bin_size {wildcards.bin_size} \
                                                        --bw_folder {params.bw_folder} \
                                                        --overall_phase_amp_table {input.overall_phase_amp} \
                                                        --out_table_amp {output.amplitude} \
                                                        --out_table_phase {output.phase} \
                                                        --out_table_mu {output.mu} \
                                                        --chr {wildcards.chr} \
                                                        --sigma {wildcards.sigma}
        """

rule plot_time_space_fourier_transform:
    input:
        amplitude=os.path.join(config['outfold_space_time_fourier'],"Bins_{bin_size}bp_sigma_{sigma}Mb","{chr}_amp_pos_k.csv"),
        phase=os.path.join(config['outfold_space_time_fourier'],"Bins_{bin_size}bp_sigma_{sigma}Mb","{chr}_phase_pos_k.csv"),
        mu=os.path.join(config['outfold_space_time_fourier'],"Bins_{bin_size}bp_sigma_{sigma}Mb","{chr}_mu_pos.csv")
    output:
        fig=os.path.join(config['outfold_fig'],"space_time_fourier","Bins_{bin_size}bp_sigma_{sigma}Mb","{chr}_amp_pos_k.pdf")
    shell:
        """
        python scripts/plot_time_space_fourier_transform.py --bin_size {wildcards.bin_size} \
                                                            --sigma {wildcards.sigma} \
                                                            --table_amp {input.amplitude} \
                                                            --table_phase {input.phase} \
                                                            --table_mu {input.mu} \
                                                            --outfig {output.fig} \
                                                            --chr {wildcards.chr}       
        """

##--------------------------------------##
##  Kalman Filter                       ##
##--------------------------------------##

rule noise_model:
    input:
        in_tables=expand(os.path.join(config['outfold_binned'],"NormCoverage_3p_bin{{bin_size}}bp_{chr}.csv"),chr=config['Chromosomes']['GRCm39'])
    output:
        noise_model_parameters=os.path.join(config['outfold_binned'],"Noise_model_parameters_{bin_size}bp.csv"),
        fig=os.path.join(config['outfold_fig'],"noise_model_{bin_size}bp.pdf")
    shell:
        """
        python scripts/noise_model.py --bin_size {wildcards.bin_size} \
                                      --out_table {output.noise_model_parameters} \
                                      --out_fig {output.fig} \
                                      --in_tables {input.in_tables}
        """

rule Kalman_gene_smoothing:
    input:
        tables=expand(os.path.join(config['outfold_binned'],"{sample}","NormCoverage_3p_{strand}_bin{{bin_size}}bp.bw"), sample=config['Samples'], strand=config['Strands']), # only for dependency
        gtf=os.path.join(config['outfold_genome'],config['Species'],"gene_protein_coding.gtf"),
        noise_model_parameters=os.path.join(config['outfold_binned'],"Noise_model_parameters_{bin_size}bp.csv")
    output:
        os.path.join(config['outfold_kalman'],"Gene_Q_1e-3","Kalman_Smoothing_{bin_size}bp.hdf5")
    params:
        bw_folder=config['outfold_binned']
    threads:
        26
    shell:
        """
        python scripts/Kalman_filter_on_genes.py --bin_size {wildcards.bin_size} \
                                                 --bw_folder {params.bw_folder} \
                                                 --gtf {input.gtf} \
                                                 --noise_model_parameters {input.noise_model_parameters} \
                                                 --out_hdf5 {output} \
                                                 --n_threads {threads}
        """

rule plot_Kalman_gene_smoothing:
    input:
        in_hdf5=os.path.join(config['outfold_kalman'],"Gene_Q_1e-3","Kalman_Smoothing_{bin_size}bp.hdf5"),
        gtf=os.path.join(config['outfold_genome'],config['Species'],"gene_protein_coding.gtf")
    output:
        fig=os.path.join(config['outfold_fig'],"kalman","Gene_Q_1e-3","scatter_Kmax_R2_LL_strand_{bin_size}bp.pdf")
    params:
        bw_folder=config['outfold_binned']
    shell:
        """
        python scripts/plot_kalman_on_genes.py --bin_size {wildcards.bin_size} \
                                               --bw_folder {params.bw_folder} \
                                               --gtf {input.gtf} \
                                               --in_hdf5 {input.in_hdf5} \
                                               --out_fig {output.fig}
        """

rule Kalman_smoothing:
    input:
        expand(os.path.join(config['outfold_binned'],"{sample}","NormCoverage_3p_{strand}_bin{{bin_size}}bp.bw"), sample=config['Samples'], strand=config['Strands']) # only for dependency
    output:
        ll=os.path.join(config['outfold_kalman'],"ll_{win_size}_{d_win}_{chr}_{bin_size}bp.csv"),
        mu=os.path.join(config['outfold_kalman'],"mu_{win_size}_{d_win}_{chr}_{bin_size}bp.npy"),
        sigma=os.path.join(config['outfold_kalman'],"sigma_{win_size}_{d_win}_{chr}_{bin_size}bp.npy"),
        r2=os.path.join(config['outfold_kalman'],"r2_{win_size}_{d_win}_{chr}_{bin_size}bp.npy")
    params:
        bw_folder=config['outfold_binned']
    threads:
        12
    shell:
        """
        python scripts/Kalman_filter_on_for_waves.py --bin_size {wildcards.bin_size} \
                                                     --bw_folder {params.bw_folder} \
                                                     --chr {wildcards.chr} \
                                                     --win_size {wildcards.win_size} \
                                                     --d_win {wildcards.d_win} \
                                                     --out_ll {output.ll} \
                                                     --out_mu {output.mu} \
                                                     --out_sigma {output.sigma} \
                                                     --out_r2 {output.r2}
        """
        

rule plot_Kalman_smoothing:
    input:
        ll=os.path.join(config['outfold_kalman'],"ll_{win_size}_{d_win}_{chr}_{bin_size}bp.csv"),
        mu=os.path.join(config['outfold_kalman'],"mu_{win_size}_{d_win}_{chr}_{bin_size}bp.npy"),
        sigma=os.path.join(config['outfold_kalman'],"sigma_{win_size}_{d_win}_{chr}_{bin_size}bp.npy"),
        r2=os.path.join(config['outfold_kalman'],"r2_{win_size}_{d_win}_{chr}_{bin_size}bp.npy")
    output:
        fig=os.path.join(config['outfold_fig'],"kalman","ll_{win_size}_{d_win}_{chr}_{bin_size}bp.pdf")
    params:
        bw_folder=config['outfold_binned']
    shell:
        """
        python scripts/plot_kalman_filter.py --bin_size {wildcards.bin_size} \
                                             --bw_folder {params.bw_folder} \
                                             --chr {wildcards.chr} \
                                             --win_size {wildcards.win_size} \
                                             --d_win {wildcards.d_win} \
                                             --in_ll {input.ll} \
                                             --in_mu {input.mu} \
                                             --in_sigma {input.sigma} \
                                             --in_r2 {input.r2} \
                                             --out_ll {output.fig}
        """

##--------------------------------------##
##  Track Hub                           ##
##--------------------------------------##

