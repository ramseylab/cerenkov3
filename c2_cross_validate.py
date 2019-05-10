import pandas as pd
import numpy as np
from xgboost import XGBClassifier
from sklearn.model_selection import cross_validate
from locus_sampling.cross_validation import BalancedGroupKFold
from locus_sampling.scoring import avg_rank_scorer2
from util_report import save_x_validate


if __name__ == "__main__":
    _RANDOM_STATE = 1337

    snp_feat_path = "./cerenkov3_data/vertex/SNP/osu18_cerenkov_feat_mat_plus_group_size.tsv"
    snp_group_path = "./cerenkov3_data/vertex/SNP/osu18_groups.tsv"

    feat_df = pd.read_csv(snp_feat_path, sep="\t")
    group_map = pd.read_csv(snp_group_path, sep="\t", usecols=["name", "group_id"])
    
    Xyg = feat_df.merge(group_map, how="left", on="name")
    y = Xyg.loc[:, "label"]
    g = Xyg.loc[:, "group_id"]
    X = Xyg.drop(columns=["label", "group_id", "name"])

    # class_balance = len(y) / sum(y) - 1  # n_negative / n_positive
    rare_event_rate = sum(y) / len(y)

    cv = BalancedGroupKFold(n_splits=5, slop_allowed=0.5, random_state=_RANDOM_STATE)

    xgb_clf = XGBClassifier(verbosity=0, objective='binary:logistic', booster='gbtree', n_jobs=-1, max_delta_step=0, scale_pos_weight=1, base_score=rare_event_rate, random_state=_RANDOM_STATE)

    scoring = {'AUPRC': 'average_precision', 
               'AUROC': 'roc_auc', 
               'AVGRANK': avg_rank_scorer2(groups=g)}
    score_names = scoring.keys()
    
    """
    g_xgb_hp_grid <- g_make_hyperparameter_grid_list(list(eta=c(0.1),
                                                        nrounds=c(40),
                                                        gamma=c(10),
                                                        lambda=c(1),
                                                        subsample=c(1),
                                                        base_score=g_class_count_frac_positive,
                                                        scale_pos_weight=c(1.0),
                                                        max_depth=c(7)))
    """
    param_dist = dict(max_depth=7,
                      learning_rate=0.1,
                      n_estimators=40, 
                      gamma=10,
                      scale_pos_weight=1,
                      base_score=rare_event_rate,
                      subsample=1)

    xgb_clf = xgb_clf.set_params(**param_dist)

    n_repeats = 10

    def _run():
        for i in range(0, n_repeats):
            cv = BalancedGroupKFold(n_splits=5, slop_allowed=0.5, random_state=_RANDOM_STATE + i)
            result_dict = cross_validate(estimator=xgb_clf, X=X, y=y, groups=g, scoring=scoring, cv=cv, return_train_score=True)
            yield result_dict
        
    result_dict_list = list(_run())

    output_dir = "./experiment_result"
    output_tag = "c2_cross_validate"
    save_x_validate(result_dict_list, output_dir, output_tag, score_names=score_names, train_score=True)