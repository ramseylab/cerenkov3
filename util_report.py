import os
import pickle
import pandas as pd
from itertools import product


def parse_hp_search_result(result_dict, score_names, train_score=True, n_splits=None, compact=False):
    # Step 1: Params
    param_df = pd.DataFrame(result_dict['params'])

    # Step 2: Ranks
    # ranks have keys like "rank_test_roc_auc"
    # No rank on training scores
    rank_keys = ["rank_test_{}".format(score_name) for score_name in score_names]
    rank_dict = {key: value for key, value in result_dict.items() if key in rank_keys}
    rank_df = pd.DataFrame(data=rank_dict)

    # Step 3: Score Stats
    # score stats have keys like "mean_test_roc_auc", with a general form "{data_name}_{stat_name}_{score_name}"
    data_names = ["test", "train"] if train_score else ["test"]
    stat_names = ["std", "mean"]
    
    score_stat_keys = ["_".join(triple) for triple in (product(stat_names, data_names, score_names))]
    score_stat_dict = {key: value for key, value in result_dict.items() if key in score_stat_keys}
    score_stat_df = pd.DataFrame(data=score_stat_dict)

    # Step 4: Scores per split
    if n_splits is not None:
        # split scores have keys like "split0_train_roc_auc"
        split_names = ["split{}".format(i) for i in range(n_splits)]
        
        split_score_keys = ["_".join(triple) for triple in (product(split_names, data_names, score_names))]
        split_score_dict = {key: value for key, value in result_dict.items() if key in split_score_keys}
        split_score_df = pd.DataFrame(data=split_score_dict)

    # Step 5: Combine
    if n_splits is not None:
        report_df = pd.concat([param_df, rank_df, score_stat_df, split_score_df], axis=1)
    else:
        report_df = pd.concat([param_df, rank_df, score_stat_df], axis=1)
    
    # sort by the first rank column
    report_df.sort_values(by=rank_keys[0], axis=0, ascending=True, inplace=True)
    
    if compact:
        # Remove columns with only one unique value
        report_df = report_df.loc[:, report_df.apply(pd.Series.nunique) != 1]
    
    return report_df

# def save_hp_search(search, folder, tag, score_names, train_score=True, n_splits=None, compact=False):
#     search_pickle_fn = '{}_search.pickle'.format(tag)
#     report_tsv_fn = '{}_report.tsv'.format(tag)

#     # Saving the "search" object
#     with open(os.path.join(folder, search_pickle_fn), 'wb') as handle:
#         pickle.dump(search, handle)

#     report_df = parse_hp_search_result(search.cv_results_, score_names, train_score, n_splits, compact)
#     report_df.to_csv(os.path.join(folder, report_tsv_fn), sep="\t", index=False)

def save_hp_search(search_list, folder, tag, score_names, train_score=True, n_splits=None, compact=False):
    object_pickle_fn = '{}_hps_object.pickle'.format(tag)
    report_tsv_fn = '{}_hps_report.tsv'.format(tag)

    # Saving the "search_list" object
    with open(os.path.join(folder, object_pickle_fn), 'wb') as handle:
        pickle.dump(search_list, handle)

    result_dict_list = [search.cv_results_ for search in search_list]

    def _parse():
        nonlocal result_dict_list

        for index, result_dict in enumerate(result_dict_list):
            report_df = parse_hp_search_result(result_dict, score_names, train_score, n_splits, compact)
            # .insert is an inplace operation by default
            report_df.insert(loc=0, column="repeat", value=index)

            yield report_df

    report_df_list = list(_parse())
    if len(report_df_list) == 1:
        report_df = report_df_list[0]
    else:
        report_df = pd.concat(report_df_list, axis=0)
    
    report_df.to_csv(os.path.join(folder, report_tsv_fn), sep="\t", index=False)

def save_x_validate(result_dict_list, folder, tag, score_names, train_score=True):
    result_dict_fn = '{}_xv_result.pickle'.format(tag)
    report_tsv_fn = '{}_xv_report.tsv'.format(tag)

    ## Saving the "result_dict" object
    with open(os.path.join(folder, result_dict_fn), 'wb') as handle:
        pickle.dump(result_dict_list, handle)

    def _parse():
        nonlocal result_dict_list

        data_names = ["test", "train"] if train_score else ["test"]    
        stat_keys = ["_".join(pair) for pair in (product(data_names, score_names))]  # E.g. "train_average_precision"

        for index, result_dict in enumerate(result_dict_list):
            sub_result_dict = dict((key, result_dict[key]) for key in stat_keys if key in result_dict)
            report_df = pd.DataFrame(data=sub_result_dict)
            # .insert is an inplace operation by default
            report_df.insert(loc=0, column="repeat", value=index)

            yield report_df

    report_df_list = list(_parse())
    if len(report_df_list) == 1:
        report_df = report_df_list[0]
    else:
        report_df = pd.concat(report_df_list, axis=0)
    
    report_df.to_csv(os.path.join(folder, report_tsv_fn), sep="\t", index=False)
    