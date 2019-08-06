import pandas as pd
import numpy as np
from xgboost import XGBClassifier
from sklearn.model_selection import RandomizedSearchCV
from locus_sampling.cross_validation import BalancedGroupKFold
from locus_sampling.scoring import avg_rank_scorer2
from util_report import save_hp_search
from c3_experiment_skeletion import prepare_for_experiment, _RANDOM_STATE

def get_param_dist():
    sgn_weight_dist = dict(w_4DGp=[0.0, 0.1, 0.3, 1.0, 3.0],
                           w_4DGt=[0.0, 0.1, 0.3, 1.0, 3.0],
                           w_GTEx=[0.0, 0.1, 0.3, 1.0, 3.0],
                           w_TFBS=[0.0, 0.1, 0.3, 1.0, 3.0],
                           w_NG=[0.0, 0.1, 0.3, 1.0, 3.0],
                           w_coexp=[0.0, 0.1, 0.3, 1.0, 3.0],
                           w_hn=[0.0, 0.1, 0.3, 1.0, 3.0],
                           w_bg=[0.0, 0.1, 0.3, 1.0, 3.0])
    n2v_param_dist = dict(d=[1, 2, 3, 4],
                          l=[2, 3, 4, 5, 6],
                          r=[3, 6, 9, 12, 15],
                          k=[2, 3, 4],
                          p=[0.125, 0.25, 0.5, 1, 2, 4, 8],
                          q=[0.125, 0.25, 0.5, 1, 2, 4, 8],
                          w=[True])
    clf_param_dist = dict(classifier__max_depth=[2, 3, 4, 5, 6, 7, 8, 9, 10],
                          classifier__learning_rate=[0.1],
                          classifier__n_estimators=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100], 
                          classifier__gamma=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                          classifier__subsample=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
                          classifier__colsample_bytree=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])

    # sgn_weight_dist = dict(w_4DGp=[3.0],
    #                        w_4DGt=[1.0],
    #                        w_GTEx=[3.0],
    #                        w_TFBS=[1.0],
    #                        w_NG=[0.3],
    #                        w_coexp=[0.1],
    #                        w_hn=[0.0],
    #                        w_bg=[0.3])
    # n2v_param_dist = dict(d=[4,5],
    #                       l=[6],
    #                       r=[6],
    #                       k=[3],
    #                       p=[4],
    #                       q=[0.125],
    #                       w=[True])
    # clf_param_dist = dict(classifier__max_depth=[8],
    #                       classifier__learning_rate=[0.1],
    #                       classifier__n_estimators=[80],
    #                       classifier__gamma=[9],
    #                       classifier__subsample=[1.0],
    #                       classifier__colsample_bytree=[0.5])

    param_dist = {**sgn_weight_dist, **n2v_param_dist, **clf_param_dist}

    return param_dist


if __name__ == "__main__":
    snp_feat_path = "./cerenkov3_data/vertex/SNP/osu19_cerenkov_feat_mat_plus_group_size.tsv"
    snp_id_path = "./cerenkov3_classifier/INT_ID_EDGELIST/SNP_INT_ID.tsv"
    snp_group_path = "./cerenkov3_data/vertex/SNP/osu19_groups.tsv"

    X_INT_ID, y, g, c3c = prepare_for_experiment(snp_feat_path, snp_id_path, snp_group_path)

    param_dist = get_param_dist()
    
    scoring = {'AUPRC': 'average_precision', 
               'AUROC': 'roc_auc', 
               'AVGRANK': avg_rank_scorer2(groups=g)}
    score_names = scoring.keys()

    n_iter_search = 240
    n_repeats = 1
    n_splits = 5

    def _run():
        for i in range(0, n_repeats):
            cv = BalancedGroupKFold(n_splits=n_splits, slop_allowed=0.5, random_state=_RANDOM_STATE + i)
            search = RandomizedSearchCV(c3c, param_distributions=param_dist,
                                        n_iter=n_iter_search, cv=cv, scoring=scoring, random_state=_RANDOM_STATE, refit=False, return_train_score=True)
            search.fit(X_INT_ID, y, groups=g)
            yield search

    search_list = list(_run())

    output_dir = "./experiment_result"
    output_tag = "c3_hyperparamter_search"
    save_hp_search(search_list, output_dir, output_tag, score_names=score_names, train_score=True, n_splits=n_splits, compact=False)