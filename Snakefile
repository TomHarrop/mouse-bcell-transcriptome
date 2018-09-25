#!/usr/bin/env python3

import pandas

###########
# GLOBALS #
###########

run_info_file = 'data/SraRunTable.txt'

sra_container = 'shub://TomHarrop/singularity-containers:sra_2.9.2'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
biopython_container = 'shub://TomHarrop/singularity-containers:biopython_1.72'
hmmer_container = 'shub://TomHarrop/singularity-containers:hmmer_3.2.1'
porechop_container = 'shub://TomHarrop/singularity-containers:porechop_0.2.3'
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.1'
clustalo_container = 'shub://TomHarrop/singularity-containers:clustalo_1.2.4'
bioc_container = 'shub://TomHarrop/singularity-containers:bioconductor_3.7'

########
# MAIN #
########

run_info = pandas.read_csv(run_info_file)

minion_runs = run_info[run_info['Model'] == 'MinION']

# split the run info into name_to_url dict
col_to_sn = minion_runs.to_dict()['SampleName']
col_to_url = minion_runs.to_dict()['download_path']

name_to_url = {}
for key in col_to_sn:
    name_to_url[col_to_sn[key]] = col_to_url[key]

all_samples = sorted(set(name_to_url.keys()))

#########
# RULES #
#########

rule target:
    input:
        'output/alignments/p450_all_trim2.faa'

rule trimal1:
    input:
        'output/alignments/p450_all_trim_align2.faa'
    output:
        'output/alignments/p450_all_trim2.faa'
    log:
        'output/logs/trimal2.log'
    threads:
        1
    singularity:
        clustalo_container
    shell:
        'trimal '
        '-in {input} '
        '-out {output} '
        '-strictplus '
        '&> {log}'

# align2
rule clustalo2:
    input:
        fasta = 'output/alignments/p450_all_trim1.faa',
        hmm = 'data/p450.hmm'
    output:
        'output/alignments/p450_all_trim_align2.faa'
    log:
        'output/logs/clustalo2.log'
    threads:
        50
    singularity:
        clustalo_container
    shell:
        'clustalo '
        '--threads={threads} '
        '--in {input.fasta} '
#        '--hmm-in {input.hmm} '
        '--dealign '
        '--out {output} '
        '&> {log}'


# trim1
rule trimal:
    input:
        'output/alignments/p450_all.faa'
    output:
        'output/alignments/p450_all_trim1.faa'
    log:
        'output/logs/trimal.log'
    threads:
        1
    singularity:
        clustalo_container
    shell:
        'trimal '
        '-in {input} '
        '-out {output} '
        '-gappyout '
        '&> {log}'

# align
rule clustalo:
    input:
        fasta = expand('output/p450_transcripts/rf{rf}.fasta',
                       rf=[1, 2, 3, 4, 5, 6]),
        mm = 'data/mm_p450s.fasta',
        hmm = 'data/p450.hmm'
    output:
        'output/alignments/p450_all.faa'
    log:
        'output/logs/clustalo.log'
    threads:
        50
    singularity:
        clustalo_container
    shell:
        'cat {input.fasta} {input.mm} | '
        'clustalo '
        '--threads={threads} '
        '--in - '
        '--hmm-in {input.hmm} '
        '--dealign '
        '--out {output} '
        '&> {log}'


# run bbmap to extract protein sequences
rule filter_by_name:
    input:
        p450_reads = 'output/hmmer/P450_read_names_rf{rf}.txt',
        fasta = 'output/translated/rf{rf}.fasta'
    output:
        'output/p450_transcripts/rf{rf}.fasta'
    threads:
        1
    log:
        'output/logs/filter_by_name/rf{rf}.log'
    singularity:
        bbduk_container
    shell:
        'filterbyname.sh '
        'ignorejunk '
        'in={input.fasta} '
        'out={output} '
        'include=t '
        'names={input.p450_reads} '
        '&> {log}'

# run R to filter reads
rule filter_p450_reads:
    input:
        hmmer_results = expand('output/hmmer/rf{rf}.txt',
                               rf=[1, 2, 3, 4, 5, 6])
    output:
        expand('output/hmmer/P450_read_names_rf{rf}.txt',
               rf=[1, 2, 3, 4, 5, 6])
    params:
        out_prefix = 'output/hmmer/P450_read_names_rf'
    log:
        'output/logs/filter_p450_reads.log'
    singularity:
        r_container
    script:
        'src/filter_p450_reads.R'


rule hmmer:
    input:
        fasta = 'output/translated/rf{rf}.fasta',
        hmm = 'data/p450.hmm'
    output:
        tbl = 'output/hmmer/rf{rf}.txt'
    threads:
        50
    log:
        'output/logs/rf{rf}_hmmer.log'
    singularity:
        hmmer_container
    shell:
        'hmmsearch '
        '--cpu {threads} '
        '--tblout {output.tbl} '
        '{input.hmm} '
        '{input.fasta} '
        '&> {log}'

rule translate:
    input:
        fasta = 'output/fasta/bcell_reads_filtered.fasta'
    output:
        rf1 = 'output/translated/rf1.fasta',
        rf2 = 'output/translated/rf2.fasta',
        rf3 = 'output/translated/rf3.fasta',
        rf4 = 'output/translated/rf4.fasta',
        rf5 = 'output/translated/rf5.fasta',
        rf6 = 'output/translated/rf6.fasta'
    threads:
        1
    singularity:
        biopython_container
    script:
        'src/translate_np_reads.py'

rule reformat:
    input:
        expand('output/porechop/{sample_name}.fastq',
               sample_name=all_samples)
    output:
        fa = 'output/fasta/bcell_reads_filtered.fasta',
        bhist = 'output/fasta/bcell_reads_filtered_bhist.txt',
        qhist = 'output/fasta/bcell_reads_filtered_qhist.txt',
    log:
        'output/logs/reformat.log'
    singularity:
        bbduk_container
    shell:
        'cat {input} | '
        'reformat.sh '
        'in=stdin.fastq '
        'int=f '
        'qin=33 '
        'out={output.fa} '
        'bhist={output.bhist} '
        'qhist={output.qhist} '
        'minlength=1000 '
        'maxlength=4000 '
        'qtrim=r '
        'trimq=7 '
        '2> {log}'

rule porechop:
    input:
        'output/fastq/{sample_name}/{sample_name}.fastq'
    output:
        'output/porechop/{sample_name}.fastq'
    threads:
        2
    log:
        'output/logs/porechop/{sample_name}.log'
    singularity:
        porechop_container
    shell:
        'porechop '
        '-i {input} '
        '-o {output} '
        '-t {threads} '
        '&> {log}'

rule dump_fastq:
    input:
        'output/SRAs/{sample_name}.sra'
    output:
        fq = 'output/fastq/{sample_name}/{sample_name}.fastq',
        tmpdir = temp(directory('output/fastq/tmp_{sample_name}'))
    threads:
        48
    params:
        outdir = 'output/fastq/{sample_name}'
    log:
        'output/logs/dump_fastq/{sample_name}.log'
    singularity:
        sra_container
    shell:
        'fasterq-dump '
        '--outfile {wildcards.sample_name}.fastq '
        '--outdir {params.outdir} '
        '--temp {output.tmpdir} '
        '--threads {threads} '
        '--details '
        '--split-files '
        '--log-level 5 '
        '{input} '
        '&> {log} '

rule download_sra:
    output:
        temp('output/SRAs/{sample_name}.sra')
    params:
        url = lambda wildcards: name_to_url[wildcards.sample_name]
    threads:
        1
    log:
        'output/logs/download_sra/{sample_name}.log'
    shell:
        'wget '
        '-O {output} '
        '{params.url} '
        '&> {log}'
