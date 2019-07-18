import pandas as pd
from scipy.stats import ttest_rel

# for col in test_cols:
#     print(col)
#     print("============")
#     print("C1 + LS (mean): {}; C3 (mean): {}".format(c2_perf[col].mean(), c3_perf[col].mean()))
#     print("t-statistic: {}; p-value: {}".format(*ttest_rel(c2_perf[col], c3_perf[col])))
#     print("\n")

def t_test_report(perf_df_a, tag_a, perf_df_b, tag_b, metric_cols):
    for col in metric_cols:
        report = dict(A=tag_a, B=tag_b, metric=col, 
                      mean_A=perf_df_a[col].mean(),
                      std_A=perf_df_a[col].std(),
                      mean_B=perf_df_b[col].mean(),
                      std_B=perf_df_b[col].std())
        
        t, p = ttest_rel(perf_df_a[col], perf_df_b[col])
        report["t-statistic"] = t
        report["p-value"] = p

        yield report

if __name__ == "__main__":
    metric_cols = ["test_AUPRC", "test_AUROC", "test_AVGRANK"]

    c1_perf = pd.read_csv("./experiment_result/c1_cross_validate_xv_report.tsv", sep="\t", usecols=metric_cols)
    c1_n2v_perf = pd.read_csv("./experiment_result/c3_minus_LS_cross_validate_xv_report.tsv", sep="\t", usecols=metric_cols)
    c1_ls_perf = pd.read_csv("./experiment_result/c2_cross_validate_xv_report.tsv", sep="\t", usecols=metric_cols)
    # C3 = C1 + LS + N2V
    c1_ls_n2v_perf = pd.read_csv("./experiment_result/c3_cross_validate_xv_report.tsv", sep="\t", usecols=metric_cols)

    report_list = list(t_test_report(c1_perf, "C1", c1_n2v_perf, "C1+N2V", metric_cols))
    report_df_1 = pd.DataFrame(report_list)

    report_list = list(t_test_report(c1_perf, "C1", c1_ls_perf, "C1+LS", metric_cols))
    report_df_2 = pd.DataFrame(report_list)

    report_list = list(t_test_report(c1_ls_perf, "C1+LS", c1_ls_n2v_perf, "C1+LS+N2V", metric_cols))
    report_df_3 = pd.DataFrame(report_list)

    report_df = pd.concat([report_df_1, report_df_2, report_df_3], axis=0)
    report_df.to_csv("t_test.tsv", sep="\t", index=False)