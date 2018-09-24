#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO

# GLOBALS

read_file = snakemake.input['fasta']
rf1_file = snakemake.output['rf1']
rf2_file = snakemake.output['rf2']
rf3_file = snakemake.output['rf3']
rf4_file = snakemake.output['rf4']
rf5_file = snakemake.output['rf5']
rf6_file = snakemake.output['rf6']

# MAIN

# stream and translate, a bit slower than reading once but saves memory
SeqIO.write((x.translate(id=x.id)
             for x in SeqIO.parse(read_file, 'fasta')),
            rf1_file,
            'fasta')

SeqIO.write((x[1:].translate(id=x.id)
             for x in SeqIO.parse(read_file, 'fasta')),
            rf2_file,
            'fasta')

SeqIO.write((x[2:].translate(id=x.id)
             for x in SeqIO.parse(read_file, 'fasta')),
            rf3_file,
            'fasta')

SeqIO.write((x.reverse_complement().translate(id=x.id)
             for x in SeqIO.parse(read_file, 'fasta')),
            rf4_file,
            'fasta')

SeqIO.write((x.reverse_complement()[1:].translate(id=x.id)
             for x in SeqIO.parse(read_file, 'fasta')),
            rf5_file,
            'fasta')

SeqIO.write((x.reverse_complement()[2:].translate(id=x.id)
             for x in SeqIO.parse(read_file, 'fasta')),
            rf6_file,
            'fasta')
