import sys
import os
import argparse
import subprocess
import math
from random import choice, randint, sample
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import namedtuple

Pop = namedtuple("Pop", [ "pop", "alt_freq", "snp_n" ])


def read_hetbed(het_bed):
    pops = []
    with open(het_bed, 'r') as fin:
        for line in fin.readlines():
            pop_id, alt_freq, snp_n = line.strip("\n").split(" ")
            pops.append(Pop(pop_id, float(alt_freq), int(snp_n)))

    return pops


def gen_vcf(pops, fasta):

    nucs = ['A', 'T', 'C', 'G']

    recs = list(SeqIO.parse(fasta, "fasta"))

    fout = open("cparv_vars.vcf", 'w')
    fout.write("##fileformat=VCFv4.2\n")
    fout.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency among genotypes, for each ALT allele, in the same order as listed">\n')
    fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    pos_box = []
    pop_bins = [[pop, near_split(pop.snp_n, len(recs))] for pop in pops]

    for i, r in enumerate(recs):
        seq = str(r.seq)
        rec_range = list(range(1, len(seq)+1))

        for pop, bins in pop_bins:
            bin = bins[i]

            for b in range(bin):
                p = choice(rec_range)
                ref = seq[p-1]
                alt = nucs[(nucs.index(ref)+1)%4]
                print(f"{r.id}\t{p}\t.\t{ref}\t{alt}\t.\t.\tAF={pop.alt_freq}", file=fout)
                rec_range.remove(p)

    fout.close()

    return "./cparv_vars.vcf"


def make_fastq(args, vcf):

    jeter_line = f"java -jar {args.jvarkit} biostar404363 -o {args.prefix}.jeter.sam -p {vcf} {args.bam}"
    sam2fq_line = f"samtools fastq -1 {args.prefix}_1.fq -2 {args.prefix}_2.fq -n {args.prefix}.jeter.sam"
    gzip_line = f"gzip {args.prefix}_1.fq {args.prefix}_2.fq"

    print(jeter_line)
    subprocess.call(jeter_line, shell=True)
    print(sam2fq_line)
    subprocess.call(sam2fq_line, shell=True)
    print(gzip_line)
    subprocess.call(gzip_line, shell=True)


def near_split(x, num_bins):

    quotient, remainder = divmod(x, num_bins)
    return [quotient + 1] * remainder + [quotient] * (num_bins - remainder)


def parseArgs(argv):
    """ simple argument parser
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('fasta', action='store', help='A sequence file in FASTA format to introduce variants into and simulate reads from.')
    parser.add_argument('jvarkit', action='store', help='Path to jvarkit.jar executable.')
    parser.add_argument('popBed', action='store', help='BED file containing definitions for each population to build this dataset.')
    parser.add_argument('bam', action='store', help='BAM file to introduce SNPs into.')

    parser.add_argument('-p', '--prefix', action='store', default="varSim", help='Prefix for output files. Default=varSim.')

    args = parser.parse_args(argv)
    return args


def main(argv):
    """ Main function
    """
    args = parseArgs(argv)

    pops = read_hetbed(args.popBed)

    vcf = gen_vcf(pops, args.fasta)

    make_fastq(args, vcf)


if __name__=="__main__":
    main(sys.argv)
