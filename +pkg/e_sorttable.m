function [T]=e_sorttable(T)
T = sortrows(T,'p_val_adj','ascend');
T = sortrows(T,'avg_log2FC','ascend');