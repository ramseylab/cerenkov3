import pandas as pd
from xgboost import XGBClassifier
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV, cross_validate
from cerenkov3_classifier.cerenkov3_classifier import Cerenkov3Classifier
from cerenkov3_classifier.node2vec_manager import N2vManager
from cerenkov3_classifier.extra_feat_manager import ExtraFeatManager
from cerenkov3_classifier.network_manager import SgnManager
from locus_sampling.cross_validation import BalancedGroupKFold

_RANDOM_STATE = 1337


def prepare_for_experiment(fixed_cv=False):
    snp_feat_path = "./cerenkov3_data/vertex/SNP/osu18_cerenkov_feat_mat_plus_group_size.tsv"
    snp_id_path = "./cerenkov3_classifier/INT_ID_EDGELIST/SNP_INT_ID.tsv"
    snp_group_path = "./cerenkov3_data/vertex/SNP/osu18_groups.tsv"

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
        #                   base_score=[rare_event_rate/100],
        #                 #   min_child_weight=[3/rare_event_rate],
        #                   subsample=[1])
    base_clf = XGBClassifier(verbosity=0, objective='binary:logistic', booster='gbtree', n_jobs=-1, scale_pos_weight=1, base_score=rare_event_rate, random_state=_RANDOM_STATE)

    c3c = Cerenkov3Classifier(sgn_manager, n2v_manager, ext_manager, base_clf)

    if fixed_cv:
        partition_table_path = "./cerenkov3_data/vertex/SNP/osu18_replications_fold_assignments.tsv"
        partition_table = pd.read_csv(partition_table_path, sep="\t")
        partition_table = Xyg.loc[:, ["name"]].merge(partition_table, how="inner", on="name")  # align partition_table to X

        return X_INT_ID, y, g, c3c, partition_table
    else:
        return X_INT_ID, y, g, c3c, None


def x_validate(sgn_weight_dist, n2v_param_dist, clf_param_dist, scoring, n_repeats=1, n_splits=5, fixed_cv=False):    
    X_INT_ID, y, g, c3c, partition_table = prepare_for_experiment(fixed_cv=fixed_cv)

    param_dist = {**sgn_weight_dist, **n2v_param_dist, **clf_param_dist}
    c3c = c3c.set_params(**param_dist)

    def _run():
        for i in range(0, n_repeats):
            if fixed_cv:
                cv = FixedReplicatedKFold(n_splits=n_splits, partition_table=partition_table, repli_colname="replication{}".format(i+1))
            else:
                cv = BalancedGroupKFold(n_splits=n_splits, slop_allowed=0.5, random_state=_RANDOM_STATE + i)
            
            result_dict = cross_validate(estimator=c3c, X=X_INT_ID, y=y, groups=g, scoring=scoring, cv=cv, return_train_score=True)
            yield result_dict
        
    result_dict_list = list(_run())

    return result_dict_list

def random_hp_search(sgn_weight_dist, n2v_param_dist, clf_param_dist, scoring, n_iter_search=None, n_repeats=1, n_splits=5, fixed_cv=False)
    X_INT_ID, y, g, c3c, partition_table = prepare_for_experiment(fixed_cv=fixed_cv)

    param_dist = {**sgn_weight_dist, **n2v_param_dist, **clf_param_dist}

    def _run():
        for i in range(0, n_repeats):
            if fixed_cv:
                cv = FixedReplicatedKFold(n_splits=n_splits, partition_table=partition_table, repli_colname="replication{}".format(i+1))
            else:
                cv = BalancedGroupKFold(n_splits=n_splits, slop_allowed=0.5, random_state=_RANDOM_STATE + i)

            search = RandomizedSearchCV(c3c, param_distributions=param_dist,
                                        n_iter=n_iter_search, cv=cv, scoring=scoring, random_state=_RANDOM_STATE, refit=False, return_train_score=True)
            search.fit(X_INT_ID, y, groups=g)
            yield search

    search_list = list(_run())

    return search_list

    
def grid_hp_search(sgn_weight_dist, n2v_param_dist, clf_param_dist, scoring, n_repeats=1, n_splits=5, fixed_cv=False):
    X_INT_ID, y, g, c3c, partition_table = prepare_for_experiment(fixed_cv=fixed_cv)
    
    param_dist = {**sgn_weight_dist, **n2v_param_dist, **clf_param_dist}
    
    def _run():
        for i in range(0, n_repeats):
            if fixed_cv:
                cv = FixedReplicatedKFold(n_splits=n_splits, partition_table=partition_table, repli_colname="replication{}".format(i+1))
            else:
                cv = BalancedGroupKFold(n_splits=n_splits, slop_allowed=0.5, random_state=_RANDOM_STATE + i)

            search = GridSearchCV(c3c, param_grid=param_dist, cv=cv, scoring=scoring, refit=False, return_train_score=True)
            search.fit(X_INT_ID, y, groups=g)
            yield search
    
    search_list = list(_run())

    return search_list
