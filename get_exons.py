
import argparse
import gffutils
import pyfaidx
import sqlite3
import os



def main():

    os.system("rm exons.fasta")


    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--transcript_id', required=True)
    args = parser.parse_args()


    #try:
        #db = gffutils.create_db('../cross_regulation_pipeline/data/input_data/gencode.v19.annotation.gtf', 'gtf.db', disable_infer_genes=True, disable_infer_transcripts=True)
    #except sqlite3.OperationalError:
        #db = gffutils.FeatureDB('gtf.db')
    #fasta = pyfaidx.Fasta('GRCh37.primary_assembly.genome.fa')

    #for cds in db.features_of_type('cds', order_by='start'):
        #print(cds.sequence(fasta))


    with open(args.gtf, "r") as gtf:
        gtf_data = gtf.readlines()

    exons = set()

    for line in gtf_data:
        if not line.startswith('#'):
            chrom, _, feature, start, stop, _, strand, _, anno = line.strip().split("\t")
            if feature == 'CDS':
                gene_id = anno.strip().split('; ')[0].split(' ')[1].strip('"')
                transcript_id = anno.strip().split('; ')[1].split(' ')[1].strip('"').split(".")[0]
                gene_type = anno.strip().split('; ')[2].split(' ')[1].strip('"')
                gene_status = anno.strip().split('; ')[3].split(' ')[1].strip('"')
                gene_name = anno.strip().split('; ')[4].split(' ')[1].strip('"')
                transcript_type = anno.strip().split('; ')[5].split(' ')[1].strip('"')
                transcript_name = anno.strip().split('; ')[7].split(' ')[1].strip('"')
                exon_number = anno.strip().split('; ')[8].split(' ')[1].strip('"')
                exon_name = anno.strip().split('; ')[9].split(' ')[1].strip('"')

                if transcript_type == 'protein_coding':
                    with open('exons.bed','a') as bed_file:
                        if strand=='+':
                        #bed_file.write(chrom, start, stop, exon_name, exon_number, strand, sep="\t")
                            bed_file.writelines("\t".join([chrom, start, str(int(stop)+1), str(gene_name+'_'+transcript_id+'_'+exon_name+'_'+exon_number),exon_number, strand]))
                            bed_file.writelines('\n')
                        elif strand=='-':
                            bed_file.writelines("\t".join([chrom, str(int(start)-1), stop, str(gene_name+'_'+transcript_id+'_'+exon_name+'_'+exon_number),exon_number, strand]))
                            bed_file.writelines('\n')
                    bed_file.close()
#exons.add(tuple([transcript_name, exon_name, exon_number, chrom, start, stop, strand]))
    os.system("sort -V -k1,1 -k2,2 exons.bed > sorted_exons.bed")
    os.system("bedtools getfasta -fi {0} -bed sorted_exons.bed  -name -s -fo exons.fasta".format(args.fasta))

    os.system("rm exons.bed")
    os.system("rm sorted_exons.bed")

    genes_exons = set()
    with open('exons.fasta','r') as exons:
        for line1,line2 in zip(exons, exons):
            if line1.startswith('>'):
                gene=line1.split('_')[0].strip('>')
                transcript = line1.split('_')[1]
                exon_name = line1.split('_')[2]
                exon_number = line1.split('_')[3].split(":")[0]
                chrom = line1.split('_')[3].split(":")[2]
                coords = line1.split('_')[3].split(":")[3].split("(")[0]
                strand = line1.split('_')[3].split(":")[3].split("(")[1].split(")")[0]

                seq = line2.strip('\n')
                genes_exons.add(tuple([gene, transcript, exon_name, exon_number, chrom, coords, strand,seq]))


    transcript_dict = {}
    for item in genes_exons:
        if item[1] == args.transcript_id:
            if item[1] in transcript_dict:
                transcript_dict[item[2]].append([item[3],item[7]])
            else:
                transcript_dict[item[2]] = [item[3],item[7]]

    sorted_transcript_dict = {k: v for k, v in sorted(transcript_dict.items(), key=lambda item:int(item[1][0]))}
    print(sorted_transcript_dict)

    seq = ""
    for key, value in sorted_transcript_dict.items():
        seq += value[1]
        seq += '|'
    print(args.transcript_id)
    print(seq)

    table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

    protein =""

    i=0
    while len(seq)>=i+3:
        codon = seq[i:i + 3]
        if len(codon) == 3:
            if codon[0]=='|':
                protein += '|'
                i = i+1
                codon = seq[i:i + 3]
                protein += table[codon]
            elif codon[1]=='|':
                protein+='|'
                codon=seq[i]+seq[i+2:i+4]
                protein += table[codon]
                i=i+1
            elif codon[2] == '|':
                if len(seq)>i+3:
                    codon = seq[i:i+2]+seq[i+3]
                    protein+=table[codon]
                    protein+='|'
                    i=i+1
                else:
                    break
            else:
                protein += table[codon]
        else:
            continue
        i=i+3
    print(protein)
if __name__ == '__main__':
    main()
