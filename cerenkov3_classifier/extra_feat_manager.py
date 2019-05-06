import pandas as pd

class ExtraFeatManager(object):
    def __init__(self, feat_df, feat_id_colname, feat_label_colname, id_map, orig_id_colname, int_id_colname):
        if feat_id_colname not in feat_df:
            raise ValueError("Column '{}' does not exist in feat_df".format(feat_id_colname))

        if feat_label_colname not in feat_df:
            raise ValueError("Column '{}' does not exist in feat_df".format(feat_label_colname))

        if int_id_colname in feat_df:
            raise ValueError("Column '{}' cannot exist in feat_df. Please rename in id_map.".format(int_id_colname))

        if len(id_map.columns) != 2:
            raise ValueError("Dataframe id_map should contain exactly 2 columns. Got {}.".format(len(id_map)))

        if orig_id_colname not in id_map:
            raise ValueError("Column '{}' does not exist in id_map".format(orig_id_colname))

        if int_id_colname not in id_map:
            raise ValueError("Column '{}' does not exist in id_map".format(int_id_colname))

        self.feat_df = feat_df
        self.feat_id_colname = feat_id_colname
        self.feat_label_colname = feat_label_colname
        
        self.id_map = id_map
        self.orig_id_colname = orig_id_colname
        self.int_id_colname = int_id_colname

    def load_feat(self):
        """
        This function guarantees to return a dataframe with columns `feat_id_colname`, `feat_label_colname` and `int_id_colname`
        """
        if self.feat_id_colname == self.orig_id_colname:
            return self.feat_df.merge(self.id_map, how="left", on=self.feat_id_colname)
        else:
            _feat_df = self.feat_df.merge(self.id_map, how="left", left_on=self.feat_id_colname, right_on=self.orig_id_colname)
            _feat_df.drop(columns=self.orig_id_colname, inplace=True)
            return _feat_df
