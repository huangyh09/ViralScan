# This program is to make GTF file for pseudo-gene based on the (viral) genome.
# It will treat the whole genome as one gene, for the purpose to get STAR-solo
# running.

# Author: Yuanhua Huang
# Date: 2020-01-01

import os
import sys
from optparse import OptionParser


class LightFasta:
    """This class is to load and to handle fasta file, which
    must be small size, e.g. less than 100M in plain text.
    This because it load the full data; big fasta file will
    take a lot of memenory and a long time to load. For big
    fasta files, you could use pysam.FastaFile. """
    def __init__(self, fasta_file):
        fid = open(fasta_file,"r")
        all_lines = fid.readlines()
        fid.close()
        seq, self.ref, self.seq = "", [], []
        for line in all_lines:
            line = line.split("\n")[0]
            if len(self.ref) == 0 and line[0] != ">": continue
            if line.startswith(">"):
                self.ref.append(line[1:])
                if seq == "": continue
                self.seq.append(seq)
                seq = ""
            else:
                seq = seq + line
        self.seq.append(seq)


def fasta_write(fid, seq, ref, length=60):
    """write sequence into a fasta file
    
    Parameter:
    ----------
    fid: a file object, the file for writing to
    seq: a string, the sequence for writing
    ref: a string, the reference id of the sequence
    length: an int, the length of sequence each line
    Return:
    -------
    None
    """
    fid.writelines(">" + ref + "\n")
    i = -1
    for i in range(int(len(seq) / length)):
        fid.writelines(seq[i*length:(i+1)*length] + "\n")
    if (i+1)*length < len(seq):
        fid.writelines(seq[(i+1)*length:] + "\n")
    return None
    
    
def main():
    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--fasta_file", "-f", dest="fasta_file",
        help="The input fasta file for viral genomes.")
    parser.add_option("--out_dir", "-o", dest="out_dir", default=None, 
        help="The directory for output files. [Default: Same as input file]")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to Viral GTF maker!")
        print("use -h or --help for help on argument.")
        sys.exit(1)
        
    if options.out_dir is None:
        out_dir = os.path.dirname(options.fasta_file)
    else:
        out_dir = options.out_dir
    
    file_basename = os.path.splitext(os.path.basename(options.fasta_file))[0]
    print("Viral GTF Maker:", file_basename)

    FASTA_FILE = LightFasta(options.fasta_file)
    FASTA_FILE.ref = ['virus' + str(x+1) + ' ' + 
                      FASTA_FILE.ref[x] for x in range(len(FASTA_FILE.ref))]
    
    ## Generate GTF
    fid = open(out_dir + "/" + file_basename + ".numbered.gtf", "w")
    for i in range(len(FASTA_FILE.ref)):
        _ref_val = FASTA_FILE.ref[i].split(" ")
        virus_num, virus_name = _ref_val[0], " ".join(_ref_val[1:])

        _line1 = [virus_num, 'virusGenome', 'gene', '1', 
                  str(len(FASTA_FILE.seq[i])), '.', '+', '.', 
                  'gene_id "%s"; gene_biotype "whole_genome"; gene_name "%s"' 
                  %(virus_num, virus_name)]
        _line2 = [virus_num, 'virusGenome', 'transcript', '1', 
                  str(len(FASTA_FILE.seq[i])), '.', '+', '.', 
                  'gene_id "%s"; transcript_id "%s"; gene_biotype "whole_genome"' 
                  %(virus_num, virus_num)]
        _line3 = [virus_num, 'virusGenome', 'exon', '1', 
                  str(len(FASTA_FILE.seq[i])), '.', '+', '.', 
                  'gene_id "%s"; transcript_id "%s"; gene_biotype "whole_genome"' 
                  %(virus_num, virus_num)]

        fid.writelines('\t'.join(_line1) + '\n')
        fid.writelines('\t'.join(_line2) + '\n')
        fid.writelines('\t'.join(_line3) + '\n')
    fid.close()
    
    
    ## Generate genome notes
    fid = open(out_dir + "/" + file_basename + ".id_notes.tsv", "w")
    for i in range(len(FASTA_FILE.ref)):
        _ref_val = FASTA_FILE.ref[i].split(" ")
        virus_num, virus_name = _ref_val[0], " ".join(_ref_val[1:])
        fid.writelines(virus_num + '\t' + virus_name + '\n')
    fid.close()
    
    
    ## Generate numbered FASTA
    fid = open(out_dir + "/" + file_basename + ".numbered.fasta", "w")
    for i in range(len(FASTA_FILE.seq)):
        fasta_write(fid, FASTA_FILE.seq[i], FASTA_FILE.ref[i], length=70)
    fid.close()
    
    

if __name__ == "__main__":
    main()

