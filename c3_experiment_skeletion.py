import pandas as pd
import numpy as np
from itertools import product
from xgboost import XGBClassifier
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV, cross_validate
from cerenkov3_classifier.cerenkov3_classifier import Cerenkov3Classifier
from cerenkov3_classifier.node2vec_manager import N2vManager
from cerenkov3_classifier.extra_feat_manager import ExtraFeatManager
from cerenkov3_classifier.network_manager import SgnManager
from util_report import save_hp_search, save_repeated_hp_searches

_RANDOM_STATE = 1337


def prepare_for_experiment():
    snp_feat_path = "./cerenkov3_data/vertex/SNP/osu18_cerenkov_feat_mat_plus_group_size.tsv"
    snp_id_path = "./cerenkov3_classifier/INT_ID_EDGELIST/SNP_INT_ID.tsv"
    snp_group_path = "./cerenkov3_data/vertex/SNP/osu18_groups.tsv"

    feat_df = pd.read_csv(snp_feat_path, sep="\t")

    label_map = feat_df.loc[:, ["name", "label"]]
    id_map = pd.read_csv(snp_id_path, sep="\t")
    group_map = pd.read_csv(snp_group_path, sep="\t", usecols=["name", "group_id"])
    
    Xyg = label_map.merge(group_map, how="left", on="name")
    Xyg = Xyg.merge(Xyg, how="left", left_on="name", right_on="ID")

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
        #                   base_score=[rare_event_rate/100],
        #                 #   min_child_weight=[3/rare_event_rate],
        #                   subsample=[1])
    base_clf = XGBClassifier(verbosity=0, objective='binary:logistic', booster='gbtree', n_jobs=-1, scale_pos_weight=1, base_score=rare_event_rate, random_state=_RANDOM_STATE)

    c3c = Cerenkov3Classifier(sgn_manager, n2v_manager, ext_manager, base_clf)

    return X_INT_ID, y, g, c3c


def x_validata(tag, sgn_weight_dist, n2v_param_dist, clf_param_dist, score_names=['average_precision', 'roc_auc']):    
    X_INT_ID, y, g, c3c = prepare_for_experiment()

    param_dist = {**sgn_weight_dist, **n2v_param_dist, **clf_param_dist}
    c3c = c3c.set_params(**param_dist)

    cv = BalancedGroupKFold(n_splits=5, slop_allowed=0.5, random_state=_RANDOM_STATE)
    scores = cross_validate(estimator=c3c, X=X_INT_ID, y=y, groups=g, scoring=score_names, cv=cv, return_train_score=True)

    save_x_validate(scores, tag, score_names)

def hp_search(tag, sgn_weight_dist, n2v_param_dist, clf_param_dist, search_type, score_names=['average_precision', 'roc_auc'], n_iter_search=None):
    param_dist = {**sgn_weight_dist, **n2v_param_dist, **clf_param_dist}
    cv = BalancedGroupKFold(n_splits=5, slop_allowed=0.5, random_state=_RANDOM_STATE)
    X_INT_ID, y, g, c3c = prepare_for_experiment()

    if search_type == "random":
        search = RandomizedSearchCV(c3c, param_distributions=param_dist,
                                    n_iter=n_iter_search, cv=cv, scoring=score_names, random_state=_RANDOM_STATE, refit=False, return_train_score=True)
        search.fit(X_INT_ID, y, groups=g)
    elif search_type == "grid":
        search = GridSearchCV(c3c, param_grid=param_dist, cv=cv, scoring=score_names, refit=False, return_train_score=True)
        search.fit(X_INT_ID, y, groups=g)
    else:
        raise ValueError("Invalid search_type. Got '{}'".format(search_type))

    tag = "{}_{}".format(tag, search_type)
    save_hp_search(search, tag, score_names)
    
def repeated_hp_searches(tag, sgn_weight_dist, n2v_param_dist, clf_param_dist, search_type, score_names=['average_precision', 'roc_auc'], n_iter_search=None, n_repeats=10):
    param_dist = {**sgn_weight_dist, **n2v_param_dist, **clf_param_dist}
    X_INT_ID, y, g, c3c = prepare_for_experiment()

    def _run():
        if search_type == "random":
            for i in range(0, n_repeats):
                cv = BalancedGroupKFold(n_splits=5, slop_allowed=0.5, random_state=_RANDOM_STATE + i)
                search = RandomizedSearchCV(c3c, param_distributions=param_dist,
                                            n_iter=n_iter_search, cv=cv, scoring=score_names, random_state=_RANDOM_STATE, refit=False, return_train_score=True)
                search.fit(X_INT_ID, y, groups=g)
                yield search
        elif search_type == "grid":
            for i in range(0, n_repeats):
                cv = BalancedGroupKFold(n_splits=5, slop_allowed=0.5, random_state=_RANDOM_STATE + i)
                search = GridSearchCV(c3c, param_grid=param_dist, cv=cv, scoring=score_names, refit=False, return_train_score=True)
                search.fit(X_INT_ID, y, groups=g)
                yield search
        else:
            raise ValueError("Invalid search_type. Got '{}'".format(search_type))

    searches = list(_run())

    tag = "{}_{}_{}".format(tag, n_repeats, search_type)
    save_repeated_hp_searches(searches, tag, score_names)

if __name__ == "__main__":
    pass