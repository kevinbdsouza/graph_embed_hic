import os
from scipy.io import loadmat


class Params:
    def __init__(self):
        self.hic_path = "/data2/hic_lstm/data/"
        self.downstream_dir = "/data2/hic_lstm/downstream"
        self.sizes_file = "chr_cum_sizes.npy"
        self.start_end_file = "starts.npy"

        self.input_file = self.hic_path + "GM12878/GM12878.hic"
        self.target_file = self.input_file
        self.label_file = None
        self.juicer_tools_path = "/data2/hic_lstm/softwares/juicer_tools.jar"
        self.dump_dir = self.hic_path + "graph"
        self.genome_size_file = "/Users/kevindsouza/Documents/UBC/PhD/Research/nucleosome/graph_embed_hic/chromosome_sizes/hg19.chrom.sizes"

        self.chr = None
        self.mode = None
        self.exp = None
        self.res = 10000
        self.order = "1"
        self.alpha = 0.5
        self.samples = 25
        self.clusters = 5
