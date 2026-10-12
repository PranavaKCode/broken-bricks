"""Final check 4 (part 3): tidy MAGMA output of the pooled-neuron gene-property run (magma_pooled.sbatch).

Reads results/final_checks/magma_pooled/<tag>/{pooled,neuron_given_microglia,microglia_given_neuron}.gsa.out
Writes results/final_checks/h3_pooled_neuron_property.csv

MAGMA reports covariate tests two-sided. As in amendment A10 the one-sided positive P is derived: P/2 if beta > 0, else 1 - P/2.
Conditional files: MODEL 1 rows are the conditioned variable first and then the other tested covariate adjusted for it.
"""
import os
import numpy as np, pandas as pd

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "final_checks", "magma_pooled")


def read(path):
    hdr = None; rows = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if hdr is None:
            hdr = p; continue
        rows.append(p[:len(hdr)])
    df = pd.DataFrame(rows, columns=hdr)
    for c in ["NGENES", "BETA", "BETA_STD", "SE", "P"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


out = []
for tag in ["primary", "neff", "noapoe"]:
    d = read(f"{R}/{tag}/pooled.gsa.out"); d["analysis"] = "marginal"; d["conditioned_on"] = ""
    out.append(d.assign(tag=tag, model=0))
    d = read(f"{R}/{tag}/neuron_given_microglia.gsa.out"); d = d[d["MODEL"].astype(str) == "1"]
    d = d[d["VARIABLE"] == "N7_Neuron_all"]; d["analysis"] = "conditional"; d["conditioned_on"] = "N7_Microglia"
    out.append(d.assign(tag=tag))
    d = read(f"{R}/{tag}/microglia_given_neuron.gsa.out"); d = d[d["MODEL"].astype(str) == "1"]
    d = d[d["VARIABLE"] == "N7_Microglia"]; d["analysis"] = "conditional"; d["conditioned_on"] = "N7_Neuron_all"
    out.append(d.assign(tag=tag))
res = pd.concat(out, ignore_index=True)
res["P_two_sided"] = res["P"]
res["P_one_sided_positive"] = np.where(res["BETA"] > 0, res["P_two_sided"] / 2, 1 - res["P_two_sided"] / 2)
res = res[["tag", "analysis", "conditioned_on", "VARIABLE", "NGENES", "BETA", "BETA_STD", "SE", "P_two_sided", "P_one_sided_positive"]]
res.to_csv(os.path.join(R, "..", "h3_pooled_neuron_property.csv"), index=False)
pd.set_option("display.width", 220)
print(res.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
