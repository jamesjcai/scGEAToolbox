import numpy as np
from typing import Optional, Dict, List

DEFAULT_POS = ["CD44", "LY6C", "KLRG1", "CTLA", "ICOS", "LAG3"]
DEFAULT_NEG = ["IL2", "TNF"]


def _get_test_data(n_genes=1000, n_samples=100, target_pos=None, random_state=42):
    random_state_seed = random_state
    random_state = np.random.default_rng(random_state)
    X = random_state.negative_binomial(20, 0.98, size=(n_genes, n_samples))
    if target_pos is None:
        target_pos = DEFAULT_POS
    gene_list = [f"pseudo_G{i}" for i in range(n_genes - len(target_pos))] + target_pos
    return {"random_state": random_state_seed,
            "X": X, "gene_list": gene_list}


def _get_assigned_bins(data_avg: np.ndarray,
                       cluster_len: int,
                       n_bins: int) -> np.ndarray:
    assigned_bin = np.zeros(shape=(cluster_len, ), dtype=np.int32)  # (G,)
    bin_size = cluster_len / n_bins
    for i_bin in range(n_bins):
        assigned_bin[(assigned_bin == 0) &
                     (data_avg <= data_avg[int(np.round(bin_size * i_bin))])] = i_bin
    return assigned_bin


def _get_ctrl_use(assigned_bin: np.ndarray,
                  gene_arr,
                  target_dict,
                  n_ctrl,
                  random_state) -> List[str]:
    selected_bins = list(set(assigned_bin[np.in1d(gene_arr, target_dict["Pos"])]))
    genes_in_same_bin = gene_arr[np.in1d(assigned_bin, selected_bins)]
    ctrl_use = list()
    for _ in range(len(target_dict["Pos"])):
        ctrl_use.extend(random_state.choice(genes_in_same_bin, n_ctrl))
    return list(set(ctrl_use))


def cell_cycle_score(X,
                     gene_list: List[str],
                     target_dict: Optional[Dict[str, List[str]]] = None,
                     n_bins: int = 25,
                     n_ctrl: int = 5,
                     random_state: int = 42):
    random_state = np.random.default_rng(random_state)
    if target_dict is None:
        target_dict = {"Pos": DEFAULT_POS,
                       "Neg": DEFAULT_NEG}
    else:
        target_dict = {k: [i.upper() for i in v] for k, v in target_dict.items()}

    if len(set(gene_list) & set(target_dict["Pos"])) == 0:
        raise ValueError('No feature genes found in gene_list.')

    gene_list = [i.upper() for i in gene_list]
    X = np.log((X / np.nansum(X, axis=0)) * 1e4 + 1)

    cluster_len = X.shape[0]
    data_avg = X.mean(axis=1)
    sort_arg = np.argsort(data_avg)
    data_avg = data_avg[sort_arg]
    gene_list = np.array(gene_list)[sort_arg]
    X = X[sort_arg, :]

    assigned_bin = _get_assigned_bins(data_avg, cluster_len, n_bins)
    used_ctrl = _get_ctrl_use(assigned_bin, gene_list, target_dict,
                              n_ctrl, random_state)
    ctrl_score = X[np.in1d(gene_list, used_ctrl), :].mean(axis=0).T
    features_score = X[np.in1d(gene_list, used_ctrl), :].mean(axis=0).T
    return features_score - ctrl_score


if __name__ == '__main__':
    test_data = _get_test_data()
    print(cell_cycle_score(**test_data))