import pandas as pd

class ExtraFeatManager(object):
    def __init__(self, feat_path, snp_id_path, id_colname, int_id_colname, label_colname):
        self.feat_path = feat_path
        self.id_map_path = id_map_path

        self.id_colname = id_colname
        self.int_id_colname = int_id_colname
        self.label_colname = label_colname

    def load_feat(self):
        # `id_colname` must exist in this feature matrix
        snp_feat_df = pd.read_csv(self.feat_path, sep="\t")

        snp_id_map = pd.read_csv(self.id_map_path, sep="\t")
        snp_id_map.columns = [self.id_colname, self.int_id_colname]

        return snp_feat_df.merge(snp_id_map, how="left", on=self.id_colname)