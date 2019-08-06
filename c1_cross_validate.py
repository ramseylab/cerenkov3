import pandas as pd
import numpy as np
from xgboost import XGBClassifier
from sklearn.model_selection import cross_validate
from locus_sampling.cross_validation import FixedReplicatedKFold
from locus_sampling.scoring import avg_rank_scorer2
from util_report import save_x_validate


if __name__ == "__main__":
    _RANDOM_STATE = 1337

    snp_feat_path = "./cerenkov3_data/vertex/SNP/osu19_cerenkov_feat_mat.tsv"
    snp_group_path = "./cerenkov3_data/vertex/SNP/osu19_groups.tsv"

    feat_df = pd.read_csv(snp_feat_path, sep="\t")
    group_map = pd.read_csv(snp_group_path, sep="\t", usecols=["name", "group_id"])
    
    Xyg = feat_df.merge(group_map, how="left", on="name")
    y = Xyg.loc[:, "label"]
    g = Xyg.loc[:, "group_id"]
    X = Xyg.drop(columns=["label", "group_id", "name"])

    # class_balance = len(y) / sum(y) - 1  # n_negative / n_positive
    rare_event_rate = sum(y) / len(y)

    xgb_clf = XGBClassifier(verbosity=0, objective='binary:logistic', booster='gbtree', n_jobs=-1, max_delta_step=0, scale_pos_weight=1, 
                            base_score=rare_event_rate, random_state=_RANDOM_STATE)
    param_dist = dict(max_depth=7,
                      learning_rate=0.1,
                      n_estimators=40, 
                      gamma=10,
                      scale_pos_weight=1,
                      base_score=rare_event_rate,
                      subsample=1)
    xgb_clf = xgb_clf.set_params(**param_dist)

    scoring = {'AUPRC': 'average_precision', 
               'AUROC': 'roc_auc', 
               'AVGRANK': avg_rank_scorer2(groups=g)}
    score_names = scoring.keys()
    
    n_repeats = 10
    n_splits = 5

    partition_table = pd.read_csv("./cerenkov3_data/vertex/SNP/osu18_replications_fold_assignments_wo_chr5_30.tsv", sep="\t")
    partition_table = Xyg.loc[:, ["name"]].merge(partition_table, how="inner", on="name")  # align partition_table to X

    def _run():
        for i in range(0, n_repeats):
            repli_colname = "replication{}".format(i+1)
            cv = FixedReplicatedKFold(n_splits=n_splits, partition_table=partition_table, repli_colname=repli_colname)
            result_dict = cross_validate(estimator=xgb_clf, X=X, y=y, groups=g, scoring=scoring, cv=cv, return_train_score=True)
            yield result_dict

    # def _run():
    #     for i in range(0, n_repeats):
    #         cv = BalancedGroupKFold(n_splits=n_splits, slop_allowed=0.5, random_state=_RANDOM_STATE + i)
    #         result_dict = cross_validate(estimator=xgb_clf, X=X, y=y, groups=g, scoring=scoring, cv=cv, return_train_score=True)
    #         yield result_dict
        
    result_dict_list = list(_run())

    output_dir = "./experiment_result"
    output_tag = "c1_cross_validate"
    save_x_validate(result_dict_list, output_dir, output_tag, score_names=score_names, train_score=True)