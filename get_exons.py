
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
            if feature == 'exon':
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
                        #bed_file.write(chrom, start, stop, exon_name, exon_number, strand, sep="\t")
                        bed_file.writelines("\t".join([chrom, start, stop, str(gene_name+'_'+transcript_id+'_'+exon_name+'_'+exon_number),exon_number, strand]))
                        bed_file.writelines('\n')
                    bed_file.close()
#exons.add(tuple([transcript_name, exon_name, exon_number, chrom, start, stop, strand]))
    os.system("sort -V -k1,1 -k2,2 exons.bed > sorted_exons.bed")
    os.system("bedtools getfasta -fi GRCh37.primary_assembly.genome.fa -bed sorted_exons.bed  -name -s -fo exons.fasta")

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

    sequence = ""
    for key, value in sorted_transcript_dict.items():
        sequence += value[1]
        sequence += '|'
    print(args.transcript_id)
    print(sequence)


if __name__ == '__main__':
    main()
