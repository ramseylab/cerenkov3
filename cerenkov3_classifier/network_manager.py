import os
import time
import json
from collections import OrderedDict
import pandas as pd
from sklearn.model_selection import ParameterGrid
import multiprocessing
import logging

logging.basicConfig(filename=__name__ + ".log",
                    filemode='a',   # append
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)


class SgnManager:
    """
    Sgn == SNP-Gene Network
    """
    def __init__(self, edgelist_dir, network_dir, network_fn_prefix="Sgn", multi_edge_policy=sum):
        self.edgelist_dir = edgelist_dir
        self.network_dir = network_dir
        self.network_fn_prefix = network_fn_prefix
        self.multi_edge_policy = multi_edge_policy

        self.DEFAULT_EDGELIST_FN_MAP = {"fn_4DGp": "SNP_x_4DGenome_promoter.INT_ID.edgelist",
                                        "fn_4DGt": "SNP_x_4DGenome_TSS.INT_ID.edgelist",
                                        "fn_GTEx": "SNP_x_GTEx.INT_ID.edgelist",
                                        "fn_TFBS": "SNP_x_EncodeTFBS.INT_ID.edgelist",
                                        "fn_NG": "SNP_x_NearestGene.INT_ID.edgelist",
                                        "fn_coexp": "Coexpedia.INT_ID.edgelist",
                                        "fn_hn": "HumanNet.INT_ID.edgelist",
                                        "fn_bg": "BioGRID.INT_ID.edgelist"}

    def _load_edgelist(self, edgelist_fn, weight=None):
        edgelist_path = os.path.join(self.edgelist_dir, edgelist_fn)
        edgelist_df = pd.read_csv(edgelist_path, sep="\t", header=None, names=["node_a", "node_b"], dtype=int)

        if weight is not None:
            edgelist_df = edgelist_df.assign(weight=weight)

        return edgelist_df

    def locate_sgn(self, w_4DGp, w_4DGt, w_GTEx, w_TFBS, w_NG, w_coexp, w_hn, w_bg):
        weights = [w_4DGp, w_4DGt, w_GTEx, w_TFBS, w_NG, w_coexp, w_hn, w_bg]
        network_fn = "{prefix}_{weights}_{multi_edge_policy}.tsv".format(prefix=self.network_fn_prefix,
                                                                         weights="_".join([str(w) for w in weights]),
                                                                         multi_edge_policy=self.multi_edge_policy.__name__)
        network_path = os.path.join(self.network_dir, network_fn)

        return network_path

    def generate_sgn(self, w_4DGp, w_4DGt, w_GTEx, w_TFBS, w_NG, w_coexp, w_hn, w_bg):
        network_path = self.locate_sgn(w_4DGp, w_4DGt, w_GTEx, w_TFBS, w_NG, w_coexp, w_hn, w_bg)

        if not os.path.exists(network_path):
            edgelist_config = OrderedDict([
                (self.DEFAULT_EDGELIST_FN_MAP["fn_4DGp"], w_4DGp),
                (self.DEFAULT_EDGELIST_FN_MAP["fn_4DGt"], w_4DGt),
                (self.DEFAULT_EDGELIST_FN_MAP["fn_GTEx"], w_GTEx),
                (self.DEFAULT_EDGELIST_FN_MAP["fn_TFBS"], w_TFBS),
                (self.DEFAULT_EDGELIST_FN_MAP["fn_NG"], w_NG),
                (self.DEFAULT_EDGELIST_FN_MAP["fn_coexp"], w_coexp),
                (self.DEFAULT_EDGELIST_FN_MAP["fn_hn"], w_hn),
                (self.DEFAULT_EDGELIST_FN_MAP["fn_bg"], w_bg),
            ])

            logger.info("composing edgelists into a SG-Network...")
            logger.info("config = {}".format(json.dumps(edgelist_config)))

            edgelists = [self._load_edgelist(path, weight)
                 for path, weight in edgelist_config.items()]
            network = pd.concat(edgelists, axis=0, ignore_index=True)

            if "weight" in network:
                network = network.groupby(["node_a", "node_b"]).aggregate(
                    self.multi_edge_policy).reset_index()
            else:  # unweighted network
                network.drop_duplicates(keep='first', inplace=True)

            network.to_csv(network_path, sep="\t", header=False, index=False)

            logger.info("SG-Network saved to {}".format(network_path))
        else:
            logger.info("SG-Network already exists, no need to generate. Location: {}".format(network_path))

        return network_path

    def _generate_sgn(self, weight_dict):
        return self.generate_sgn(**weight_dict)

    def batch_generate_sgn(self, weight_grid):
        n_cells = len(ParameterGrid(weight_grid))
        n_cores = multiprocessing.cpu_count()
        if n_cells > n_cores - 1:
            logger.info("Got {} parameter cell(s); has {} core(s); starting {} process(es)".format(n_cells, n_cores, n_cores - 1))
            pool = multiprocessing.Pool(n_cores - 1)
        else:
            logger.info("Got {} parameter cell(s); has {} core(s); starting {} process(es)".format(n_cells, n_cores, n_cells))
            pool = multiprocessing.Pool(n_cells)

        t_start = time.perf_counter()

        pool.map(self._generate_sgn, ParameterGrid(weight_grid))
        pool.close()  # No more work
        pool.join()  # Wait for completion

        t_stop = time.perf_counter()

        logger.info("Elapsed time: %.1f [sec]" % (t_stop - t_start))

# if __name__ == "__main__":
#     sgnm = SgnManager(edgelist_dir="./INT_ID_EDGELIST", network_dir="./INT_ID_NETWORK")

#     weight_grid = {"w_4DGp": [1.0],
#                    "w_4DGt": [1.0],
#                    "w_GTEx": [1.0],
#                    "w_TFBS": [1.0],
#                    "w_NG": [1.0],
#                    "w_coexp": [1.0],
#                    "w_hn": [1.0],
#                    "w_bg": [1.0, 2.0]}

#     sgnm.batch_generate_sgn(weight_grid)
