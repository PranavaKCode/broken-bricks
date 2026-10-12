#!/usr/bin/env python
"""V2: merge the analyst's judgments into the audit sample and compute precision with Wilson 95% CIs.

The judgments were made by the analysis agent (an AI model, not a human domain expert) after reading each mention string in its sentence of
the title/abstract. 'correct' = the string refers to the gene or its protein/peptide product (Abeta -> APP follows PubTator3's own convention);
'incorrect' = the string refers to something else (different entity, chemical, complex, or probe); 'unclear' = family-level or drug-class
mention that cannot be pinned to this gene. The CSV is left editable so a human can re-judge.
For rows with verifier = rejected, 'true_rejection' = the mention string is absent from the title/abstract (plain case-insensitive search)."""
import os, math
from pathlib import Path
import polars as pl

RES = Path(os.path.expanduser("~/orcd/pool/mantis-films/research/lit/out/results"))
d = pl.read_csv(RES / "ner_audit_50_sample.csv")
J = {
 0: ("correct", "S6K1 = ribosomal protein S6 kinase B1"), 1: ("incorrect", "PGD2 is prostaglandin D2 (a lipid), not the PGDS gene"),
 2: ("correct", ""), 3: ("correct", "JIP2 = JNK-interacting protein 2"), 4: ("correct", ""), 5: ("correct", ""),
 6: ("unclear", "generic 'synaptotagmin' (family), assigned to SYT1"), 7: ("correct", ""), 8: ("correct", "tPA named as a thrombolytic drug; entity is PLAT"),
 9: ("incorrect", "AP4 = adaptor protein complex 4, not transcription factor AP-4 (TFAP4)"), 10: ("correct", ""), 11: ("unclear", "beta-arrestin family (ARRB1 or ARRB2)"),
 12: ("correct", ""), 13: ("correct", ""), 14: ("correct", ""), 15: ("incorrect", "ARP = aldehyde reactive probe, a chemical reagent"),
 16: ("correct", "PKC-delta"), 17: ("correct", ""), 18: ("correct", ""), 19: ("correct", "JKAP = DUSP22"), 20: ("incorrect", "TOM = translocase of the outer membrane complex, not PRPF6"),
 21: ("correct", ""), 22: ("correct", ""), 23: ("correct", ""), 24: ("correct", "PKC alpha"),
 25: ("correct", "Abeta -> APP per PubTator3 convention (peptide product)"), 26: ("correct", ""), 27: ("correct", ""), 28: ("correct", ""), 29: ("correct", ""), 30: ("correct", ""),
 31: ("correct", ""), 32: ("correct", ""), 33: ("incorrect", "'GLP-1 receptor agonists' is a drug class / receptor, not the GCG gene"), 34: ("correct", ""),
 35: ("correct", ""), 36: ("correct", ""), 37: ("correct", "tau written as the Greek letter in the abstract"), 38: ("correct", ""), 39: ("correct", "Abeta -> APP per PubTator3 convention"),
 40: ("correct", ""), 41: ("correct", ""), 42: ("correct", ""), 43: ("correct", ""), 44: ("correct", ""), 45: ("correct", ""), 46: ("correct", ""),
 47: ("unclear", "'cholinesterase inhibition' = drug class; BCHE vs ACHE not specified"), 48: ("correct", ""), 49: ("correct", "Abeta -> APP per PubTator3 convention"),
}
for i in range(50, 60):
    J[i] = ("true_rejection", "mention string absent from title/abstract (annotation came from full text)")
rows = d.to_dicts()
for i, r in enumerate(rows):
    r["judgment"], r["note"] = J[i]
out = pl.DataFrame(rows)
out.write_csv(RES / "ner_audit_50.csv")

def wilson(k, n, z=1.96):
    if n == 0: return (float("nan"), float("nan"))
    p = k / n; den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h
res = []
for name, sel in (("A_lt_100", out.filter(pl.col("stratum") == "A_lt_100")), ("A_ge_100", out.filter(pl.col("stratum") == "A_ge_100")),
                  ("all_50_accepted", out.filter(pl.col("stratum").is_in(["A_lt_100", "A_ge_100"])))):
    n = sel.height; k = int((sel["judgment"] == "correct").sum()); inc = int((sel["judgment"] == "incorrect").sum()); un = int((sel["judgment"] == "unclear").sum())
    lo, hi = wilson(k, n); lo2, hi2 = wilson(k, n - un)
    res.append(dict(stratum=name, n=n, correct=k, incorrect=inc, unclear=un, precision_unclear_as_wrong=k / n, wilson95_lo=lo, wilson95_hi=hi,
                    precision_excl_unclear=k / (n - un), wilson95_lo_excl=lo2, wilson95_hi_excl=hi2))
rj = out.filter(pl.col("stratum") == "extra_verifier_rejected"); n = rj.height; kk = int((rj["judgment"] == "true_rejection").sum())
lo, hi = wilson(n - kk, n)
res.append(dict(stratum="extra_verifier_rejected (false rejections)", n=n, correct=kk, incorrect=n - kk, unclear=0, precision_unclear_as_wrong=None, wilson95_lo=lo, wilson95_hi=hi,
                precision_excl_unclear=None, wilson95_lo_excl=None, wilson95_hi_excl=None))
pl.DataFrame(res, infer_schema_length=None).write_csv(RES / "ner_audit_summary.csv")
print(pl.DataFrame(res, infer_schema_length=None))
