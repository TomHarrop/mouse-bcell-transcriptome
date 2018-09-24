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
        expand('output/hmmer/rf{rf}.txt',
               rf=[1, 2, 3, 4, 5, 6])

# run transdecoder to get reads of interest

# run bbmap to extract filtered reads

# run R to filter reads

rule hmmer:
    input:
        fasta = 'output/translated/rf{rf}.fasta',
        hmm = 'data/p450.hmm'
    output:
        tbl = 'output/hmmer/rf{rf}.txt'
    threads:
        50
    log:
        'output/hmmer/rf{rf}_hmmer.log'
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
        expand('output/fastq/{sample_name}/{sample_name}.fastq',
               sample_name=all_samples)
    output:
        fa = 'output/fasta/bcell_reads_filtered.fasta',
        bhist = 'output/fasta/bcell_reads_filtered_bhist.txt',
        qhist = 'output/fasta/bcell_reads_filtered_qhist.txt',
    log:
        'output/fasta/reformat.log'
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

rule join_fastq:
    input:
        expand('output/fastq/{sample_name}/{sample_name}.fastq',
               sample_name=all_samples)
    output:
        'output/fastq_all/all_bcell_reads.fastq'
    shell:
        'cat {input} > {output}'

rule dump_fastq:
    input:
        'output/SRAs/{sample_name}.sra'
    output:
        fq = 'output/fastq/{sample_name}/{sample_name}.fastq',
        tmpdir = temp(directory('output/fastq/tmp_{sample_name}'))
    priority:
        1
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
