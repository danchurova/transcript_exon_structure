# transcript_exon_structure
Scripts for parsing .gtf and .fasta files to analyze transcripts' exon-intron structure

get_exons.py: get information from annotation (.gtf) and fasta (main assembly) and print sequence with "|" between exons for a transcript specified by user

Requirements: bedtools

usage: python3 get_exons.py --gtf anno.gtf --fasta data.fa --transcript_id ENST00000582730
