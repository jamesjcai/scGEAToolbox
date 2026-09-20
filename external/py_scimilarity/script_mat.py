import _scimilarity_env  # noqa: F401  (tiledb DLL search path, cp1252-safe stdout)
from scipy.sparse import csr_matrix
import h5py
import anndata
from scimilarity.utils import align_dataset, lognorm_counts
from scimilarity import CellAnnotation

f = h5py.File("X.mat", 'r')
# Raw counts, not normalised values. Normalisation belongs after the gene
# alignment below -- see the note there.
Xcounts = csr_matrix(f.get('Xcounts'))
modeldir = f['modeldir'][()]
modeldir = modeldir.tobytes().decode('utf-16')
g = f['g']
gene_names = []
for r in g:
    for ref in r:
        str_data = ''.join(chr(c[0]) for c in f[ref][:])
        gene_names.append(str_data)

target_celltypes = []
if 'tg' in f:
    g = f['tg']
    for r in g:
        for ref in r:
            str_data = ''.join(chr(c[0]) for c in f[ref][:])
            target_celltypes.append(str_data)

f.close()

adata = anndata.AnnData(X=Xcounts)
adata.obs.index = [f"cell_{i}" for i in range(adata.n_obs)]
adata.var.index = gene_names
adata.layers["counts"] = adata.X.copy()
print("Input data ready.")

model_path = modeldir
ca = CellAnnotation(model_path=model_path)
print("Model read.")

if target_celltypes != []:
    ca.safelist_celltypes(target_celltypes)

# Align first, normalise second -- the order both SCimilarity tutorials use.
# align_dataset drops every gene outside ca.gene_order and carries
# layers["counts"] through, so lognorm_counts divides by a library size taken
# over the model's own 28k gene space. Normalising before the alignment, as
# this script used to, divides by a library size that also counts genes the
# model is about to discard, which shifts every value the encoder sees by a
# data-dependent factor.
adata = align_dataset(adata, ca.gene_order)
adata = lognorm_counts(adata)
embeddings = ca.get_embeddings(adata.X)
print("Get_embeddings...done.")

# weighting=True matches the annotation tutorial: the vote among the 50
# reference neighbours is weighted by 1/distance rather than one-neighbour-
# one-vote.
predictions, nn_idxs, nn_dists, nn_stats = ca.get_predictions_knn(
    embeddings, weighting=True
)
print("Prediction...done.")
predictions.to_csv("output.csv", index=True)

# Confidence per predicted label, which the caller uses to decide what to
# trust. "vsAll" is the winning label's share of the 50 reference neighbours
# and "vs2nd" its share against the runner-up alone; the _weighted variants
# are the same counts under the 1/distance weighting actually used for the
# prediction. All of this was already computed and thrown away.
nn_stats[["vs2nd", "vsAll", "vs2nd_weighted", "vsAll_weighted",
          "min_dist", "max_dist"]].to_csv("output_stats.csv", index=True)
print("Output written.")
