"""Final check 3 (side result): nuclei per subclass in the CELLxGENE h5ad objects (used by the CNA) against the final-nuclei metadata CSV
(used for donor totals and class shares). Writes results/final_checks/cellxgene_vs_final_nuclei_counts.csv.
Needs results/cna/<label>__all_cells.parquet (committed, Day 1) and results/final_checks/donor_subclass_counts.csv."""
import os
import numpy as np, pandas as pd
from scipy import stats

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
csv = pd.read_csv(f"{R}/final_checks/donor_subclass_counts.csv")
m = csv.pivot_table(index="donor_id", columns="Subclass", values="n", aggfunc="sum", fill_value=0)
lab = {"microglia": "Immune", "astrocyte": "Astrocyte", "oligodendrocyte": "Oligodendrocyte", "opc": "OPC", "sst": "Sst", "pvalb": "Pvalb", "vip": "Vip",
       "lamp5": "Lamp5", "l23it": "L2/3 IT", "l4it": "L4 IT", "l5it": "L5 IT"}
rows = []
for k, v in lab.items():
    c = pd.read_parquet(f"{R}/cna/{k}__all_cells.parquet", columns=["donor_id"])["donor_id"].value_counts()
    a = pd.concat([c.rename("cxg"), m[v].rename("csv")], axis=1).fillna(0)
    r = a.cxg / a.csv.replace(0, np.nan)
    rows.append(dict(subclass=v, cxg=int(a.cxg.sum()), csv=int(a.csv.sum()), ratio=a.cxg.sum() / a.csv.sum(), donor_ratio_median=r.median(),
                     donor_ratio_min=r.min(), donor_ratio_max=r.max(), spearman=stats.spearmanr(a.cxg, a.csv)[0]))
d = pd.DataFrame(rows)
print(d.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
d.to_csv(f"{R}/final_checks/cellxgene_vs_final_nuclei_counts.csv", index=False)
