import pandas as pd
import numpy as np
from itertools import product
from xgboost import XGBClassifier
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
from cerenkov3_classifier.cerenkov3_classifier import Cerenkov3Classifier
from cerenkov3_classifier.node2vec_manager import N2vManager
from cerenkov3_classifier.extra_feat_manager import ExtraFeatManager
from cerenkov3_classifier.network_manager import SgnManager


def prepare_for_test():
    snp_feat_path = "../../S_OSU18_SNP/feature_matrix_osu18_2_fillna_sparse_plus_group_size.tsv"
    snp_int_id_path = "../../Septet_SNP_x_GTEx_x_4DGenome_x_NearestGene_x_TFBS_x_HumanNet_x_Coexpedia/INT_ID_EDGELIST/SNP_INT_ID.tsv"
    snp_label_group_path = "../../S_OSU18_SNP/osu18_SNP.tsv"

    snp_int_id_df = pd.read_csv(snp_int_id_path, sep="\t")
    snp_label_group_df = pd.read_csv(snp_label_group_path, sep="\t", usecols=["name", "label", "group_id"])
    
    Xyg = snp_int_id_df.merge(snp_label_group_df, how="left", on="name")
    X = Xyg.loc[:, "INT_ID"]
    y = Xyg.loc[:, "label"]
    g = Xyg.loc[:, "group_id"]

    class_balance = len(y) / sum(y) - 1  # n_negative / n_positive
    rare_event_rate = sum(y) / len(y)

    sgn_manager = SgnManager(
        edgelist_dir="../../Septet_SNP_x_GTEx_x_4DGenome_x_NearestGene_x_TFBS_x_HumanNet_x_Coexpedia/INT_ID_EDGELIST", 
        network_dir="../../Septet_SNP_x_GTEx_x_4DGenome_x_NearestGene_x_TFBS_x_HumanNet_x_Coexpedia/INT_ID_NETWORK")
    n2v_manager = N2vManager(emb_dir="../N2V_EMB", feat_dir="../N2V_FEAT", command="../../snap/examples/node2vec/node2vec",
                             int_id_path=snp_int_id_path,
                             label_path=snp_label_group_path)
    ckv_manager = CkvManager(feat_path=snp_feat_path, 
                             int_id_path=snp_int_id_path)
    base_clf = XGBClassifier(verbosity=0, objective='binary:logistic', booster='gbtree', n_jobs=6, max_delta_step=0, scale_pos_weight=1, base_score=rare_event_rate, random_state=1337,
                             min_child_weight=3 / rare_event_rate / 100)

    sncc = SgnN2vCkvClassifier(sgn_manager, n2v_manager, ckv_manager, base_clf)

    return X, y, g, sncc


def run_test(tag, sgn_weight_dist, n2v_param_dist, clf_param_dist, search_type, score_names=['average_precision', 'roc_auc'], n_iter_search=None):
    param_dist = {**sgn_weight_dist, **n2v_param_dist, **clf_param_dist}
    cv = BalancedGroupKFold(n_splits=5, slop_allowed=0.5, random_state=1337)
    X, y, g, sncc = prepare_for_test()

    if search_type == "random":
        tag = "{}_{}".format(tag, n_iter_search)
        search = RandomizedSearchCV(sncc, param_distributions=param_dist,
                                    n_iter=n_iter_search, cv=cv, scoring=score_names, random_state=1337, refit=False)
        search.fit(X, y, groups=g)
    elif search_type == "grid":
        search = GridSearchCV(sncc, param_grid=param_dist, cv=cv, scoring=score_names, refit=False)
        search.fit(X, y, groups=g)
    else:
        raise ValueError("Invalid search_type. Got '{}'".format(search_type))

    save_search(search, tag, score_names)
    
def run_repeated_tests(tag, sgn_weight_dist, n2v_param_dist, clf_param_dist, search_type, score_names=['average_precision', 'roc_auc'], n_iter_search=None, n_repeats=10):
    param_dist = {**sgn_weight_dist, **n2v_param_dist, **clf_param_dist}
    X, y, g, sncc = prepare_for_test()

    tag = "{}_{}".format(tag, n_repeats)

    def _run():
        nonlocal tag

        for i in range(0, n_repeats):
            cv = BalancedGroupKFold(n_splits=5, slop_allowed=0.5, random_state=1337+i)

            if search_type == "random":
                search = RandomizedSearchCV(sncc, param_distributions=param_dist,
                                            n_iter=n_iter_search, cv=cv, scoring=score_names, random_state=1337, refit=False, return_train_score=True)
                search.fit(X, y, groups=g)
            elif search_type == "grid":
                search = GridSearchCV(sncc, param_grid=param_dist, cv=cv, scoring=score_names, refit=False)
                search.fit(X, y, groups=g)
            else:
                raise ValueError("Invalid search_type. Got '{}'".format(search_type))

            yield search

    searches = list(_run())
    save_repeated_searches(searches, tag, score_names)

if __name__ == "__main__":
    pass