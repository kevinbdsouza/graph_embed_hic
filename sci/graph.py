import argparse
from sci.hic import HicData
from sci.Compartments import Compartments
from sci import params


def run_sci(cfg):
    myobject = HicData(cfg.res, cfg.exp, cfg.chr)
    myobject.initialize(cfg.genome_size_file)
    myobject.load_interaction_data(cfg.input_file)
    hic_graph = myobject.write_inter_chrom_graph()
    GW_metadata = myobject.get_bins_info()

    predictor = Compartments(cfg.exp, GW_metadata)
    predictor.predict_subcompartents(hic_graph,
                                     cfg.order,
                                     cfg.samples,
                                     cfg.clusters)


if __name__ == '__main__':
    cfg = params.Params()
    cfg.exp = "run1"

    for chr in range(1, 23):
        cfg.chr = chr
        run_sci(cfg)
