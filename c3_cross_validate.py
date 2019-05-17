import pandas as pd
import numpy as np
from sklearn.model_selection import cross_validate
from locus_sampling.cross_validation import FixedReplicatedKFold
from locus_sampling.scoring import avg_rank_scorer2
from util_report import save_x_validate
from c3_experiment_skeletion import prepare_for_experiment, _RANDOM_STATE

def get_param_dist():
    """
    Best parameters from c3_hyperparameter_search.py
    """
    sgn_weight_dist = dict(w_4DGp=3.0,
                           w_4DGt=0.3,
                           w_GTEx=0.3,
                           w_TFBS=0.1,
                           w_NG=0.1,
                           w_coexp=0.3,
                           w_hn=0.3,
                           w_bg=3.0)
    n2v_param_dist = dict(d=2,
                          l=6,
                          r=12,
                          k=4,
                          p=4,
                          q=8,
                          w=True)
    clf_param_dist = dict(classifier__max_depth=10,
                          classifier__learning_rate=0.1,
                          classifier__n_estimators=100,
                          classifier__gamma=10,
                          classifier__subsample=1.0,
                          classifier__colsample_bytree=0.3)

    param_dist = {**sgn_weight_dist, **n2v_param_dist, **clf_param_dist}

    return param_dist


if __name__ == "__main__":
    snp_feat_path = "./cerenkov3_data/vertex/SNP/osu18_cerenkov_feat_mat_plus_group_size.tsv"
    snp_id_path = "./cerenkov3_classifier/INT_ID_EDGELIST/SNP_INT_ID.tsv"
    snp_group_path = "./cerenkov3_data/vertex/SNP/osu18_groups.tsv"
    partition_table_path = "./cerenkov3_data/vertex/SNP/osu18_replications_fold_assignments.tsv"

    X_INT_ID, y, g, c3c, partition_table = prepare_for_experiment(snp_feat_path, snp_id_path, snp_group_path, partition_table_path)

    param_dist = get_param_dist()
    c3c = c3c.set_params(**param_dist)

    scoring = {'AUPRC': 'average_precision', 
               'AUROC': 'roc_auc', 
               'AVGRANK': avg_rank_scorer2(groups=g)}
    score_names = scoring.keys()    
    
    n_repeats = 10
    n_splits = 5

    def _run():
        for i in range(0, n_repeats):
            cv = FixedReplicatedKFold(n_splits=n_splits, partition_table=partition_table, repli_colname="replication{}".format(i+1))
            result_dict = cross_validate(estimator=c3c, X=X_INT_ID, y=y, groups=g, scoring=scoring, cv=cv, return_train_score=True)
            yield result_dict
        
    result_dict_list = list(_run())

    output_dir = "./experiment_result"
    output_tag = "c3_cross_validate"
    save_x_validate(result_dict_list, output_dir, output_tag, score_names=score_names, train_score=True)
