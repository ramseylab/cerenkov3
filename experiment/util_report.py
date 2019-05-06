import os
import pandas as pd
import pickle


def cv_report_to_df(cv_report, score_names, train_score=False, n_split=None, compact=False):
    # Step 1: Params
    param_df = pd.DataFrame(cv_report['params'])

    # Step 2: Ranks
    # ranks have keys like "rank_test_roc_auc"
    # No rank on training scores
    rank_keys = ["rank_test_{}".format(score_name) for score_name in score_names]
    rank_dict = {key: value for key, value in cv_report.items() if key in rank_keys}
    rank_df = pd.DataFrame(data=rank_dict)

    # Step 3: Score Stats
    # score stats have keys like "mean_test_roc_auc", with a general form "{data_name}_{stat_name}_{score_name}"
    data_names = ["test", "train"] if train_score else ["test"]
    stat_names = ["std", "mean"]
    
    score_stat_keys = ["_".join(triple) for triple in (product(stat_names, data_names, score_names))]
    score_stat_dict = {key: value for key, value in cv_report.items() if key in score_stat_keys}
    score_stat_df = pd.DataFrame(data=score_stat_dict)

    # Step 4: Scores per split
    if n_split is not None:
        # split scores have keys like "split0_train_roc_auc"
        split_names = ["split{}".format(i) for i in range(n_split)]
        
        split_score_keys = ["_".join(triple) for triple in (product(split_names, data_names, score_names))]
        split_score_dict = {key: value for key, value in cv_report.items() if key in split_score_keys}
        split_score_df = pd.DataFrame(data=split_score_dict)

    # Step 5: Combine
    if n_split is not None:
        report_df = pd.concat([param_df, rank_df, score_stat_df, split_score_df], axis=1)
    else:
        report_df = pd.concat([param_df, rank_df, score_stat_df], axis=1)
    
    # sort by the first rank column
    report_df.sort_values(by=rank_keys[0], axis=0, ascending=True, inplace=True)
    
    if compact:
        # Remove columns with only one unique value
        report_df = report_df.loc[:, report_df.apply(pd.Series.nunique) != 1]
    
    return report_df

def save_search(search, folder, tag, score_names, train_score, n_split, compact):
    search_pickle_fn = '{}_search.pickle'.format(tag)
    report_pickle_fn = '{}_report.pickle'.format(tag)
    report_tsv_fn = '{}_report.tsv'.format(tag)

    ## Saving the "search" object
    with open(os.path.join(folder, search_pickle_fn), 'wb') as handle:
        pickle.dump(search, handle)

    ## Saving the "report" parsed from "search.cv_results_"
    with open(os.path.join(folder, report_pickle_fn), 'wb') as handle:
        pickle.dump(search.cv_results_, handle)

    report_df = cv_report_to_df(search.cv_results_, score_names, train_score, n_split, compact)
    print(report_df)
    report_df.to_csv(os.path.join(folder, report_tsv_fn), sep="\t", index=False)

def save_repeated_searches(search_list, tag, score_names, train_score, n_split, compact):
    search_pickle_fn = '{}_search_list.pickle'.format(tag)
    report_pickle_fn = '{}_report_list.pickle'.format(tag)
    report_tsv_fn = '{}_report.tsv'.format(tag)

    ## Saving the "search_list" object
    with open(os.path.join(folder, search_pickle_fn), 'wb') as handle:
        pickle.dump(search_list, handle)

    cv_report_list = [search.cv_results_ for search in search_list]

    with open(os.path.join(folder, report_pickle_fn), 'wb') as handle:
        pickle.dump(cv_report_list, handle)

    def _parse():
        nonlocal cv_report_list

        for rep_ind, cv_report in enumerate(cv_report_list):
            report_df = cv_report_to_df(cv_report, score_names, train_score, n_split, compact)
            # .insert is an inplace operation by default
            report_df.insert(loc=0, column="repeat", value=rep_ind)

            yield report_df

    report_df_list = list(_parse())
    report_df = pd.concat(report_df_list, axis=0)
    print(report_df)
    report_df.to_csv(os.path.join(folder, report_tsv_fn), sep="\t", index=False)