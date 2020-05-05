import os
import numpy as np
from itertools import combinations
from tqdm import tqdm
import pandas as pd


class HicData:
    def __init__(self, cfg):
        self.res = cfg.res
        self.chr = cfg.chr
        self.chrom = "chr" + str(self.chr)
        self.mode = cfg.mode
        self.cfg = cfg
        self.contact_matrices = {}
        self.dExpected_counts = {}
        self.toKeep = {}
        self.GW_contact_matrix = []
        self.GW_meta_data = []
        self.GW_toKeep = []
        self.GW_location2index = {}
        self.dChrBins = {}
        self.name = cfg.name
        self.ordered_pairs = []
        self.embedding_mode = None
        self.interchrom_contact_matrix = None

    def get_chromosomes_bin_counts(self, chrsize_file):
        oF = open(chrsize_file)
        for line in oF.readlines():
            parts = line.strip().split()
            self.dChrBins[parts[0]] = int(int(parts[1]) / self.res + 1)
        return

    def initialize(self):
        self.get_chromosomes_bin_counts(self.cfg.chr_size_file)
        # generate genomewide mapping
        for i in range(1, 23):
            chrom = "chr%d" % i
            for j in range(self.dChrBins[chrom]):
                self.GW_meta_data.append((chrom, j * self.res,
                                          j * self.res + self.res))

        for i, location in enumerate(self.GW_meta_data):
            self.GW_location2index[(location[0], location[1])] = i

    def initialize_contact_matrix(self):
        # initialize cis interaction matrices

        if self.mode == "cis":
            # for chrom in self.dChrBins.keys():
            shape = self.dChrBins[self.chrom]
            self.contact_matrices[(self.chrom, self.chrom)] = np.zeros((shape, shape))

        if self.mode == "trans":
            # initialize trans matricies
            for chrom1, chrom2 in combinations(self.dChrBins.keys(), 2):
                rows = self.dChrBins[chrom1]
                cols = self.dChrBins[chrom2]
                self.contact_matrices[(chrom1, chrom2)] = np.zeros((rows, cols))
                self.ordered_pairs.append((chrom1, chrom2))

        return

    def construct_chr(self, prefix='hic', hic_res=10000):

        filepath = os.path.join(self.cfg.dump_dir, '{1}_chrm{0}_chrm{0}.txt'.format(self.chr, prefix))
        txt_data = pd.read_csv(filepath, sep='\t', header=None).values

        rows = txt_data[:, 0] / hic_res
        cols = txt_data[:, 1] / hic_res

        data = txt_data[:, 2]

        rows = rows.astype(int)
        cols = cols.astype(int)

        self.contact_matrices[(self.chrom, self.chrom)][rows, cols] = data
        self.contact_matrices[(self.chrom, self.chrom)][cols, rows] = data

        return

    def hicToContact(self, prefix='hic'):
        """ Calls juicer_tools to extract hic data into txt files """

        M = None
        if self.mode == "trans":
            for chrm1 in range(1, 3, 2):
                for chrm2 in range(2, 4, 2):
                    output_txt_path = os.path.join(self.cfg.dump_dir,
                                                   '{2}_chrm{0}_chrm{1}.txt'.format(chrm1, chrm2, prefix))

                    os.system("java -jar {0} dump observed KR {1} {2} {3} BP 100000 {4} > tmp_juicer_log".format(
                        self.cfg.juicer_tools_path, self.cfg.input_file,
                        chrm1, chrm2,
                        output_txt_path))

        else:
            output_txt_path = os.path.join(self.cfg.dump_dir, '{1}_chrm{0}_chrm{0}.txt'.format(self.chr, prefix))
            self.cfg.output_txt_path = output_txt_path
            os.system(
                "java -jar {0} dump observed KR {1} {2} {2} BP 10000 {3}".format(self.cfg.juicer_tools_path,
                                                                                 self.cfg.input_file, self.chr,
                                                                                 output_txt_path))

        self.construct_chr()

        return

    def load_interaction_data(self):

        self.initialize_contact_matrix()

        if self.mode == "trans":
            oF = open(self.cfg.input_file, "r")
            for line in tqdm(oF.readlines(), desc='Reading %s' % self.cfg.input_file):
                (chrom1, start1, end1, chrom2,
                 start2, end2, count) = line.strip().split()
                i = int(start1) / self.res
                j = int(start2) / self.res

                # Debug snnipet
                try:
                    row, col = self.contact_matrices[(chrom1, chrom2)].shape
                    if i >= row or j >= col:
                        continue
                except KeyError:
                    row, col = self.contact_matrices[(chrom2, chrom1)].shape
                    if j >= row or i >= col:
                        continue

                try:
                    self.contact_matrices[(chrom1, chrom2)][i, j] += \
                        float(count)
                except KeyError:
                    try:
                        self.contact_matrices[(chrom2, chrom1)][j, i] += \
                            float(count)
                    except KeyError:
                        print("ignoring pair: %s,"
                              "chromoses are not in"
                              "chromsomes size file") % line.strip()
        else:
            self.hicToContact()

        return

    def get_contact_matrix(self, chrom1, chrom2):
        return self.contact_matrices[(chrom1, chrom2)]

    def get_bins_info(self):
        return self.GW_meta_data

    def write_inter_chrom_graph(self):
        outfile = os.path.join(self.cfg.dump_dir, str(self.chr) + "_HiC_graph.txt")
        self.output_graph_path = outfile

        oF = open(outfile, "w")

        if self.mode == "trans":
            for (chrom1, chrom2) in tqdm(self.ordered_pairs,
                                         desc="writing HiC graph"):
                if chrom1 != chrom2:
                    if self.embedding_mode == "oe":
                        c1 = int(chrom1.strip("chr"))
                        c2 = int(chrom2.strip("chr"))

                        if c1 % 2 == 1 and c2 % 2 == 0:
                            self.dump_interchrom_block(oF, chrom1, chrom2)
                    else:
                        self.dump_interchrom_block(oF, chrom1, chrom2)
        else:
            self.dump_interchrom_block(oF, self.chrom, self.chrom)

        return outfile

    def dump_interchrom_block(self, oF, chrom1, chrom2):
        row, col = self.contact_matrices[(chrom1, chrom2)].shape
        for i in range(row):
            for j in range(col):
                node1 = self.GW_location2index[(chrom1, i * self.res)]
                node2 = self.GW_location2index[(chrom2, j * self.res)]
                value = self.contact_matrices[(chrom1, chrom2)][i, j]
                if value >= 1:
                    oF.write("%d\t%d\t%f\n" % (node1, node2, value))
                    oF.write("%d\t%d\t%f\n" % (node2, node1, value))

    def write_GW_matrix(self, mat_file, cis=True):
        rows = []
        for i in range(1, 23):
            chrom1 = "chr%d" % i
            cols = []
            for j in range(1, 23):
                chrom2 = "chr%d" % j
                if (i == j):
                    if cis:
                        cols.append(self.contact_matrices
                                    [(chrom1, chrom1)])
                    else:
                        size = self.dChrBins[chrom1]
                        cols.append(np.zeros((size, size)))
                else:
                    try:
                        cols.append(self.contact_matrices
                                    [(chrom1, chrom2)])
                    except KeyError:
                        try:
                            cols.append(self.
                                        contact_matrices
                                        [(chrom2, chrom1)]
                                        .T)
                        except KeyError:
                            print("Can not find matrix "
                                  "for%s %s, will add "
                                  "zeroes to GW matrix"
                                  % (chrom1, chrom2))
                            row = self.dChrBins[chrom1]
                            col = self.dChrBins[chrom2]
                            cols.append(np.zeros(
                                (row, col)))
            row = np.concatenate(cols, axis=1)
            rows.append(row)
        self.interchrom_contact_matrix = np.concatenate(rows, axis=0)
        np.savetxt(mat_file, self.interchrom_contact_matrix, delimiter=",")
        return
