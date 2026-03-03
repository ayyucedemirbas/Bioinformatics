import cellxgene_census
import scanpy as sc

with cellxgene_census.open_soma() as census:
    adata = cellxgene_census.get_anndata(
        census,
        organism="Homo sapiens",
        obs_value_filter="tissue == 'brain' and disease == 'normal' and cell_type == 'microglial cell'",
        column_names={"obs": ["tissue", "disease", "cell_type", "sex"]}
    )

    print(adata)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

adata.var_names = adata.var["feature_name"]
adata.var_names_make_unique()

sc.pl.umap(adata, color=['cell_type', 'AIF1'])
