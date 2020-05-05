from sci.hic import HicData
from sci.Compartments import Compartments
from sci import params
import os
import numpy as np


def run_sci(cfg):
    myobject = HicData(cfg)
    myobject.initialize()
    myobject.load_interaction_data()
    hic_graph = myobject.write_inter_chrom_graph()
    GW_metadata = myobject.get_bins_info()

    predictor = Compartments(cfg, GW_metadata)
    features = predictor.predict_subcompartents(hic_graph)

    np.save(os.path.join(cfg.dump_dir, "embedding_" + str(cfg.chr)), features)

    os.system("rm {}".format(cfg.output_txt_path))
    os.system("rm {}".format(cfg.output_graph_path))
    os.system("rm {}".format(cfg.embedding_file_path))

    return features


if __name__ == '__main__':
    cfg = params.Params()
    cfg.exp = "run1"

    for chr in range(12, 23):
        cfg.chr = chr
        run_sci(cfg)

        print("done")
