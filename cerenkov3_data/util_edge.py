import pandas as pd


def remove_duplicated_undirected_edges(edgelist_df, sort=False):
    """
    For an undirected graph, edge (v, w) is the equivalent to (w, v).
    This function removes the duplicated edges in an edgelist dataframe.
    """
    edge_set = set()

    for row in edgelist_df.to_numpy():
        edge = tuple(row)
        if (edge not in edge_set) and (tuple(reversed(row)) not in edge_set):
            edge_set.add(edge)

    new_edgelist_df = pd.DataFrame(data=edge_set, columns=edgelist_df.columns)

    if sort:
        new_edgelist_df.sort_values(by=new_edgelist_df.columns.tolist(), inplace=True)
    
    return new_edgelist_df
    