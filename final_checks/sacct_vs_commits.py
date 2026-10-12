"""Compare Slurm (sacct) job times on ORCD with the git commit times of the preregistration and amendments.

Inputs : sacct_research_jobs.csv (from `sacct -u $USER -X -P --format=JobID,JobName,Partition,Submit,Start,End,Elapsed,State`, research jobs lit-/seaad-/d1-/d2-/d3 only)
         git log of this repository
Output : results/final_checks/sacct_vs_commits.csv and a printed table. All times are EDT (ORCD and the laptop both ran on EDT).
"""
import csv, os, re, subprocess, sys
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..")
F = "%Y-%m-%dT%H:%M:%S"

def t(s):
    return datetime.strptime(s, F) if s and s not in ("None", "Unknown") else None

jobs = []
with open(os.path.join(HERE, "sacct_research_jobs.csv"), encoding="utf-8", newline="") as fh:
    for p in csv.DictReader(fh):
        jobs.append(dict(id=p["jobid"], name=p["jobname"], part=p["partition"], sub=t(p["submit"]), start=t(p["start"]), end=t(p["end"]), state=p["state"]))

log = subprocess.run(["git", "-C", REPO, "log", "--format=%h|%ad|%s", "--date=format:" + F],
                     capture_output=True, text=True).stdout.strip().splitlines()
commits = {}
for l in log:
    h, d, s = l.split("|", 2)
    commits[h] = (t(d), s)

def fam(rx):
    r = re.compile(rx)
    return sorted([j for j in jobs if r.search(j["name"]) and j["sub"] and j["sub"].year == 2026 and j["sub"].month == 10 and j["sub"].day >= 11],
                  key=lambda j: j["sub"])

def delta(a, b):
    """b - a in h:mm:ss with sign."""
    s = int((b - a).total_seconds())
    sign = "-" if s < 0 else "+"
    s = abs(s)
    return f"{sign}{s//3600}:{(s%3600)//60:02d}:{s%60:02d}"

rows = []

def check(label, commit, family_rx, expect, key="start", filt=None):
    """expect = 'after' : every job in the family must have `key` after the commit (report earliest);
       expect = 'before': every job in the family must have `key` before the commit (report latest)."""
    ct = commits[commit][0]
    js = [j for j in fam(family_rx) if (filt is None or filt(j))]
    js = [j for j in js if j[key]]
    if not js:
        rows.append([label, commit, ct, family_rx, expect, "", "", "", "no jobs"])
        return
    if expect == "after":
        j = min(js, key=lambda j: j[key])
        bad = [x for x in js if x[key] < ct]
        ok = "OK" if not bad else f"VIOLATION ({len(bad)} jobs before)"
        rows.append([label, commit, ct, family_rx, expect, f"{j['id']} {j['name']}", j[key], delta(ct, j[key]), ok])
    else:
        j = max(js, key=lambda j: j[key])
        bad = [x for x in js if x[key] > ct]
        ok = "OK" if not bad else f"VIOLATION ({len(bad)} jobs after)"
        rows.append([label, commit, ct, family_rx, expect, f"{j['id']} {j['name']}", j[key], delta(ct, j[key]), ok])

# 1 preregistration: no research job may start before the commit
check("prereg: any lit-/seaad- job starts after the commit", "ef5646b", r"^(lit|seaad)-", "after")
check("prereg: any lit-/seaad- job SUBMITTED after the commit", "ef5646b", r"^(lit|seaad)-", "after", key="sub")
# A1: logged before any metric (earliest stage A job)
check("A1: first stage-A (metric) job starts after A1", "f0da21e", r"^lit-stage[ABD]", "after")
# A2: genetics-vs-attention joins (stage B / D) must start after A2
check("A2: first stage-B job (genetics joined to attention) starts after A2", "d54de6d", r"^lit-stageB", "after")
check("A2: first stage-D job starts after A2", "d54de6d", r"^lit-stageD", "after")
check("A2: stage-A rewrite run (A_g abstract-verified) starts after A2", "d54de6d", r"^lit-stageA-q", "after")
# A2: stage-A jobs that finished before A2 (the output A2 describes)
check("A2: original stage-A run ended before A2 (the output A2 cites)", "d54de6d", r"^lit-stageA$", "before", key="end",
      filt=lambda j: j["state"] == "COMPLETED" and j["start"] is not None and j["end"] < commits["d54de6d"][0])
# A4
check("A4: first 11-subclass CNA job (d1-) submitted after A4", "9123dc3", r"^d1-", "after", key="sub")
# A6/A7
check("A6: first Day-2 job (d2-) submitted after A6/A7", "41658d0", r"^d2-", "after", key="sub")
check("A7: first Day-3 job (d3 / d3p) submitted after A6/A7", "41658d0", r"^d3", "after", key="sub")
# A8: the debug run (d2-brk-q) ends before A8; full brick jobs after
check("A8: the debug run (d2-brk-q) ended before A8", "0b7ccf0", r"^d2-brk-q$", "before", key="end")
check("A8: first full brick analysis (d2-brk-all) starts after A8", "0b7ccf0", r"^d2-brk-all$", "after")
# A8 nuance: Day-2 NAM (graph + diffusion, phenotype-free, no brick score) jobs ran once before A8 and were rerun after it
nam_rx = re.compile(r"^d2-(micro|microglia|astrocyte|oligodendrocyte|opc|sst|pvalb|vip|lamp5|l23it|l4it|l5it)-(all|nosad)$")
nam = [j for j in fam(r"^d2-") if nam_rx.match(j["name"]) and j["state"] == "COMPLETED"]
ca = commits["0b7ccf0"][0]
names = sorted({j["name"] for j in nam})
rerun_ok = [n for n in names if any(j["name"] == n and j["start"] > ca for j in nam)]
before_ct = sum(1 for j in nam if j["end"] <= ca)
rows.append(["A8 (informational): Day-2 NAM jobs, phenotype-free, no brick score: %d of %d names have a run started after A8; %d runs ended before A8" % (len(rerun_ok), len(names), before_ct),
             "0b7ccf0", ca, "d2-<subclass>-(all|nosad)", "rerun after", "n/a", ca, "n/a",
             "INFO (first NAM runs, incl. nosad, preceded A8; they hold PC scores only)"])
# A9
check("A9: first full 'all' brick analysis ended before A9", "4232120", r"^d2-brk-all$", "before", key="end",
      filt=lambda j: j["start"] and j["start"] < commits["4232120"][0])
check("A9: nosad brick analysis (d2-brk-nosad) starts after A9", "4232120", r"^d2-brk-nosad$", "after")
check("A9: techResid rerun (second d2-brk-all) starts after A9", "4232120", r"^d2-brk-all$", "after",
      filt=lambda j: j["start"] and j["start"] > commits["4232120"][0])
# results commit: last compute job ends before it
check("results e4f87cf: last Day-2/3 job ended before the results commit", "e4f87cf", r"^(d2|d3)", "before", key="end")

with open(os.path.join(HERE, "..", "results", "final_checks", "sacct_vs_commits.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["check", "commit", "commit_time_EDT", "job_family_regex", "expected", "extreme_job", "job_time", "job_minus_commit", "result"])
    w.writerows(rows)

for r in rows:
    print(f"{r[8]:<10} {r[0]:<78} commit {r[1]} {r[2]:%H:%M:%S}  job {r[5]:<34} {r[6]:%m-%d %H:%M:%S}  {r[7]}")

# bookkeeping on reruns / failures of the analysis jobs
print()
print("State counts for analysis jobs (lit-, seaad-, d1-, d2-, d3):")
from collections import Counter
c = Counter(j["state"].split()[0] for j in jobs if re.match(r"^(lit|seaad|d1|d2|d3)", j["name"]) and j["sub"].day >= 11)
print(dict(c))
print("Non-COMPLETED analysis jobs:")
for j in sorted([j for j in jobs if re.match(r"^(lit|seaad|d1|d2|d3)", j["name"]) and j["sub"].day >= 11 and not j["state"].startswith("COMPLETED")], key=lambda j: j["sub"]):
    print(f"  {j['id']} {j['name']:<24} {j['state']:<14} start {j['start']}")
print("Same-name repeats among COMPLETED jobs (name: n runs):")
cn = Counter(j["name"] for j in jobs if re.match(r"^(lit|seaad|d1|d2|d3)", j["name"]) and j["sub"].day >= 11 and j["state"] == "COMPLETED")
print({k: v for k, v in cn.items() if v > 1})
