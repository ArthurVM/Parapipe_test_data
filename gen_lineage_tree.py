import sys
import os
import argparse
import subprocess
import math
from random import choice, randint, sample
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import namedtuple, defaultdict

Pop = namedtuple("Pop", [ "pop", "alt_freq", "snp_n" ])


def gen_bed(fasta):

    nucs = ['A', 'T', 'C', 'G']

    pops = 6
    order1_step = 100000
    order2_step = 100
    l1_s = 0
    l2_s = 10000
    l3_s = 20000

    fout = open("cp_e2etree.bed", 'w')

    lindict = defaultdict(list)

    for rec in SeqIO.parse(fasta, "fasta"):
        chr = rec.id
        seq = str(rec.seq)

        for i in range(pops):
            for i_i in range(pops):
                pos = (i_i+1)*order1_step+((i)*10000)
                ref = seq[pos]
                alt = nucs[(nucs.index(ref)+1)%4]

                linline = f"SP {chr} {pos} {pos+1} {alt} L{i+1}"
                print(linline, file=fout)

                lindict[f"L{i+1}"].append(linline)

                for j in range(pops):
                    pos += order2_step
                    ref = seq[pos]
                    alt = nucs[(nucs.index(ref)+1)%4]

                    linline = f"SP {chr} {pos} {pos+1} {alt} L{i+1}.{i_i+1}"
                    print(linline, file=fout)

                    lindict[f"L{i+1}.{i_i+1}"].append(linline)

    fout.close()

    for lin, subd in lindict.items():
        with open(lin+".bed", 'w') as fout:
            print("\n".join(subd), file=fout)
            if len(lin) > 2:
                rootlin = lin.split(".")[0]
                print("\n".join(lindict[rootlin]), file=fout)


def run_ngsc():

    rootdir = "/home/amorris/BioInf/Parapipe_modtest/datasets/variant_simulation/e2e_beds/e2e_lin_defs/"
    for bed in os.listdir(rootdir):
        if bed.endswith(".bed"):
            prefix = bed.split(".bed")[0]

            if not os.path.exists(f"{prefix}_1.fq.gz"):
                runline = f"python3 /home/amorris/BioInf/ngsContrive/ngsContrive.py /home/amorris/BioInf/ngsCWD/Cp/cryptosporidium_parvum.fasta /home/amorris/BioInf/art_bin_MountRainier/ -p {prefix}_ -d 20 -s 1.0 -v {rootdir}{bed}"
                subprocess.call(runline, shell=True)
                # print(runline)


def main(fasta):
    """ Main function
    """
    # gen_bed(fasta)
    run_ngsc()


if __name__=="__main__":
    main("/home/amorris/BioInf/Parapipe_testing_modules/variant_simulation/cryptosporidium_parvum.fasta")
