import pandas as pd
from xgboost import XGBClassifier
from cerenkov3_classifier.cerenkov3_classifier import Cerenkov3Classifier
from cerenkov3_classifier.node2vec_manager import N2vManager
from cerenkov3_classifier.extra_feat_manager import ExtraFeatManager
from cerenkov3_classifier.network_manager import SgnManager
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV, cross_validate
from locus_sampling.cross_validation import BalancedGroupKFold, FixedReplicatedKFold

_RANDOM_STATE = 1337


def prepare_for_experiment(snp_feat_path, snp_id_path, snp_group_path, partition_table_path):
    # snp_feat_path = "./cerenkov3_data/vertex/SNP/osu18_cerenkov_feat_mat_plus_group_size.tsv"
    # snp_id_path = "./cerenkov3_classifier/INT_ID_EDGELIST/SNP_INT_ID.tsv"
    # snp_group_path = "./cerenkov3_data/vertex/SNP/osu18_groups.tsv"
    # partition_table_path = "./cerenkov3_data/vertex/SNP/osu18_replications_fold_assignments.tsv"

    feat_df = pd.read_csv(snp_feat_path, sep="\t")
    label_map = feat_df.loc[:, ["name", "label"]]
    id_map = pd.read_csv(snp_id_path, sep="\t")
    group_map = pd.read_csv(snp_group_path, sep="\t", usecols=["name", "group_id"])
    
    Xyg = label_map.merge(group_map, how="left", on="name")
    Xyg = Xyg.merge(id_map, how="left", left_on="name", right_on="ID")

    X_INT_ID = Xyg.loc[:, "INT_ID"]
    y = Xyg.loc[:, "label"]
    g = Xyg.loc[:, "group_id"]

    # class_balance = len(y) / sum(y) - 1  # n_negative / n_positive
    rare_event_rate = sum(y) / len(y)

    sgn_manager = SgnManager(edgelist_dir="./cerenkov3_classifier/INT_ID_EDGELIST", 
                             network_dir="./cerenkov3_classifier/INT_ID_NETWORK")
    n2v_manager = N2vManager(emb_dir="./cerenkov3_classifier/N2V_EMB", 
                             feat_dir="./cerenkov3_classifier/N2V_FEAT", 
                             command="./cerenkov3_classifier/node2vec",
                             id_map=id_map, orig_id_colname="ID", int_id_colname="INT_ID", 
                             label_map=label_map, label_id_colname="name", label_colname="label")
    ext_manager = ExtraFeatManager(feat_df=feat_df, feat_id_colname="name", feat_label_colname="label", 
                                   id_map=id_map, orig_id_colname="ID", int_id_colname="INT_ID")
    
    # The best parameters for CERENKOV2 are:
        # param_dist = dict(max_depth=[7],
        #                   learning_rate=[0.1],
        #                   n_estimators=[40], 
        #                   gamma=[10],
        #                   scale_pos_weight=[1],
        #                   base_score=[rare_event_rate],
        #                   subsample=[1])
    base_clf = XGBClassifier(verbosity=0, objective='binary:logistic', booster='gbtree', n_jobs=-1, scale_pos_weight=1, base_score=rare_event_rate, random_state=_RANDOM_STATE)

    c3c = Cerenkov3Classifier(sgn_manager, n2v_manager, ext_manager, base_clf)

    if partition_table_path:
        partition_table = pd.read_csv(partition_table_path, sep="\t")
        partition_table = Xyg.loc[:, ["name"]].merge(partition_table, how="inner", on="name")  # align partition_table to X

        return X_INT_ID, y, g, c3c, partition_table
    else:
        return X_INT_ID, y, g, c3c, None
