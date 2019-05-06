"""
/////////////////////////////////////////////////////////////////////////////

Parameters:
Input graph path (-i:)
Output graph path (-o:)
Number of dimensions. Default is 128 (-d:)
Length of walk per source. Default is 80 (-l:)
Number of walks per source. Default is 10 (-r:)
Context size for optimization. Default is 10 (-k:)
Number of epochs in SGD. Default is 1 (-e:)
Return hyperparameter. Default is 1 (-p:)
Inout hyperparameter. Default is 1 (-q:)
Verbose output. (-v)
Graph is directed. (-dr)
Graph is weighted. (-w)
Output random walks instead of embeddings. (-ow)

/////////////////////////////////////////////////////////////////////////////

Usage:
./node2vec -i:graph/karate.edgelist -o:emb/karate.emb -l:3 -d:24 -p:0.3 -dr -v
"""

import subprocess
import os
import time
import pandas as pd
import logging

logging.basicConfig(filename=__name__ + ".log",
                    filemode='a',   # append
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)


class N2vManager(object):
    def __init__(self, emb_dir, feat_dir, command, id_map, orig_id_colname, int_id_colname, label_map, label_id_colname, label_colname):
        self.emb_dir = emb_dir
        self.feat_dir = feat_dir
        self.command = command

        if len(id_map.columns) != 2:
            raise ValueError("Dataframe id_map should contain exactly 2 columns. Got {}.".format(len(id_map)))

        if orig_id_colname not in id_map:
            raise ValueError("Column '{}' does not exist in id_map".format(orig_id_colname))

        if int_id_colname not in id_map:
            raise ValueError("Column '{}' does not exist in id_map".format(int_id_colname))

        if len(label_map.columns) != 2:
            raise ValueError("Dataframe label_map should contain exactly 2 columns. Got {}.".format(len(label_map)))

        if label_id_colname not in label_map:
            raise ValueError("Column '{}' does not exist in label_map".format(label_id_colname))

        if label_colname not in label_map:
            raise ValueError("Column '{}' does not exist in label_map".format(label_colname))

        if label_colname == int_id_colname:
            raise ValueError("Param 'label_colname' cannot be the same with 'int_id_colname'. Please rename.")

        if any(colname.startswith("node2vec_") for colname in [orig_id_colname, int_id_colname, label_id_colname, label_colname]):
            raise ValueError("Any input column name cannot start with 'node2vec_'. This prefix is preserved for Node2vec features.")
        
        self.id_map = id_map
        self.orig_id_colname = orig_id_colname
        self.int_id_colname = int_id_colname

        self.label_map = label_map
        self.label_id_colname = label_id_colname
        self.label_colname = label_colname

        # `label_id_map_` will contain 3 columns: label_id_colname, label_colname, int_id_colname
        if self.label_id_colname == self.orig_id_colname:
            self.label_id_map_ = self.label_map.merge(self.id_map, how="left", on=self.label_id_colname)
        else:
            self.label_id_map_ = self.label_map.merge(self.id_map, how="left", left_on=self.label_id_colname, right_on=self.orig_id_colname)
            self.label_id_map_.drop(columns=self.orig_id_colname, inplace=True)

    def _format_command_param(self, d, l, r, k, p, q, w, i, o):
        yield "-d:{}".format(d)
        yield "-l:{}".format(l)
        yield "-r:{}".format(r)
        yield "-k:{}".format(k)
        yield "-p:{}".format(p)
        yield "-q:{}".format(q)

        if w:
            yield "-w"

        yield "-i:{}".format(i)
        yield "-o:{}".format(o)
        
    def locate_emb(self, network_fn, d, l, r, k, p, q, w):
        params = [d, l, r, k, p, q, w]

        prefix = network_fn.split(".tsv")[0]
        infix = "_".join([str(p) for p in params])
        emb_fn = "{prefix}_{infix}.emb".format(prefix=prefix, infix=infix)
        emb_path = os.path.join(self.emb_dir, emb_fn)

        return emb_path

    def generate_emb(self, network_path, d, l, r, k, p, q, w):
        network_fn = os.path.basename(network_path)
        emb_path = self.locate_emb(network_fn, d, l, r, k, p, q, w)

        if not os.path.exists(emb_path):
            logger.info("running node2vec...")
            logger.info("params: {}".format(" ".join(self._format_command_param(d=d, l=l,
                                                                                r=r, k=k, 
                                                                                p=p, q=q, w=w, 
                                                                                i=network_path, o=emb_path))))
            
            # OpenMP is automatically enabled, so all cores are used when running node2vec
            subprocess.run([self.command, *self._format_command_param(d=d, l=l,
                                                                      r=r, k=k, 
                                                                      p=p, q=q, w=w, 
                                                                      i=network_path, o=emb_path)], check=True)

            logger.info("embedding saved to {}".format(emb_path))
        else:
            logger.info("embedding already exists, no need to generate. Location: {}".format(emb_path))

        return emb_path

    def _read_emb(self, path):
        # the 1st row specifies the dimensions, "nrow ncol"
        emb_df = pd.read_csv(path, delim_whitespace=True, skiprows=1, header=None)

        col_names = [self.int_id_colname] + \
            ["node2vec_{}".format(i) for i in range(1, emb_df.shape[1])]
        emb_df.columns = col_names

        return emb_df

    def locate_feat(self, network_fn, d, l, r, k, p, q, w):
        params = [d, l, r, k, p, q, w]

        prefix = network_fn.split(".tsv")[0]
        infix = "_".join([str(p) for p in params])
        feat_fn = "{prefix}_{infix}.tsv".format(prefix=prefix, infix=infix)
        feat_path = os.path.join(self.feat_dir, feat_fn)

        return feat_path

    def generate_feat(self, network_path, d, l, r, k, p, q, w):
        network_fn = os.path.basename(network_path)
        feat_path = self.locate_feat(network_fn, d, l, r, k, p, q, w)

        if not os.path.exists(feat_path):
            logger.info("generating features from embedding...")

            emb_path = self.generate_emb(network_path, d, l, r, k, p, q, w)
            # emb_df has only INT_ID (self.int_id_colname), no rsID, no label
            emb_df = self._read_emb(emb_path)

            # now we add rsID and label back into `emb_df`
            feat_df = self.label_id_map_.merge(emb_df, how="left", on=self.int_id_colname)
            # `feat_df` guarantees to be a dataframe with columns `int_id_colname`, `label_id_colname` and `label_colname`
            feat_df.to_csv(feat_path, sep="\t", header=True, index=False)
            
            logger.info("feature saved to {}".format(feat_path))
        else:
            logger.info("feature already exists, no need to generate. Location: {}".format(feat_path))

        return feat_path

# if __name__ == "__main__":
#     id_map = pd.read_csv("./INT_ID_EDGELIST/SNP_INT_ID.tsv", sep="\t")  # columns = ["ID", "INT_ID"]
#     label_map = pd.read_csv("../cerenkov3_data/vertex/SNP/osu18_SNP_basic_info.tsv", sep="\t", usecols=["name", "label"])

#     n2vm = N2vManager(emb_dir="./N2V_EMB", feat_dir="./N2V_FEAT", command="./node2vec", 
#                       id_map=id_map, orig_id_colname="ID", int_id_colname="INT_ID",
#                       label_map=label_map, label_id_colname="name", label_colname="label")

#     network_path = "./INT_ID_NETWORK/Sgn_1.0_1.0_1.0_1.0_1.0_1.0_1.0_1.0_sum.tsv"
#     train_params = dict(
#         d=2,
#         l=2,
#         r=10,
#         k=2,
#         p=1,
#         q=1,
#         w=True
#     )

#     n2vm.generate_feat(network_path, **train_params)
