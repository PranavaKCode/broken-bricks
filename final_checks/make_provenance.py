#!/usr/bin/env python
"""Write results/provenance.csv: when the preregistration was committed, when each later commit was made, what had and had not been
looked at before the preregistration commit, and the data sources and versions.

Computed from: `git log` of this repository, final_checks/sacct_research_jobs.csv (Slurm accounting of the research jobs on ORCD, EDT),
results/data_provenance.csv (literature-study retrievals) and the public GitHub API (repository creation and first push). The SEA-AD,
GWAS, MAGMA and software rows are descriptive entries taken from the Day-0 data report; their retrieval windows are read from the sacct
table by job id.

Run from anywhere:  python final_checks/make_provenance.py      (needs network access for the GitHub API rows)
Commits made after this file was generated are not listed; regenerate to include them.
"""
import csv
import hashlib
import json
import re
import subprocess
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "results" / "provenance.csv"
EDT = timezone(timedelta(hours=-4))
UTC = timezone.utc
PREREG = "ef5646bea1ed45c144384904b3f577015140a09c"


def fmt(dt):
    """(EDT string, UTC string) for an aware datetime."""
    if dt is None:
        return "", ""
    return dt.astimezone(EDT).strftime("%Y-%m-%dT%H:%M:%S-04:00"), dt.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def edt_naive(s):
    return datetime.strptime(s, "%Y-%m-%dT%H:%M:%S").replace(tzinfo=EDT) if s and s != "None" else None


rows = []


def add(section, item, value, dt=None, source=""):
    e, u = fmt(dt)
    rows.append([section, item, value, e, u, source])


# ---------------------------------------------------------------------------------------------- git history
log = subprocess.run(["git", "-C", str(REPO), "log", "--reverse", "--format=%H%x1f%cI%x1f%s"], capture_output=True, text=True, check=True).stdout.strip().splitlines()
commits = []
for l in log:
    h, d, s = l.split("\x1f", 2)
    commits.append((h, datetime.fromisoformat(d), s))
t_prereg = next(d for h, d, s in commits if h == PREREG)


def kind(subject):
    if subject.startswith("Preregistration"):
        return "preregistration"
    if subject.startswith("Amendment") or subject.startswith("Fix timestamp in A2"):
        return "amendment"
    if subject.startswith("Prepare for public release"):
        return "release"
    if re.search(r"\b(results|Day 1:|Days 2-3|Final checks|Report inputs)\b", subject):
        return "results"
    return "code"


add("preregistration", "commit", PREREG, t_prereg, "git log; commit time is the committer time recorded in the commit")
add("preregistration", "file_unchanged_since_commit",
    "yes" if subprocess.run(["git", "-C", str(REPO), "diff", "--quiet", PREREG, "HEAD", "--", "preregistration.md"]).returncode == 0 else "no",
    None, "git diff of preregistration.md between the preregistration commit and HEAD")

proto = (REPO / "protocol" / "06-research-project.md").read_bytes()
add("preregistration", "protocol_copy_sha256", hashlib.sha256(proto).hexdigest().upper() + f" ({proto.count(b'\n')} lines)", None,
    "protocol/06-research-project.md; the same SHA-256 is recorded in the header of preregistration.md")

for h, d, s in commits:
    add(f"commit_{kind(s)}", h[:7], s, d, "git log")

head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
add("table", "commits_listed_through", head, None, "git rev-parse HEAD when this table was generated; later commits are not listed")

# ---------------------------------------------------------------------------------------------- independent clocks (GitHub)
def gh(path):
    req = urllib.request.Request("https://api.github.com/" + path, headers={"User-Agent": "broken-bricks-provenance", "Accept": "application/vnd.github+json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


repo = gh("repos/PranavaKCode/broken-bricks")
add("independent_clock", "github_repository_created", "PranavaKCode/broken-bricks", datetime.fromisoformat(repo["created_at"].replace("Z", "+00:00")),
    "GitHub API repos/PranavaKCode/broken-bricks created_at")
try:
    ev = [e for e in gh("repos/PranavaKCode/broken-bricks/events?per_page=100") if e["type"] == "CreateEvent"]
    ev.sort(key=lambda e: e["created_at"])
    add("independent_clock", "github_first_push_branch_created", ev[0]["payload"].get("ref", "main") if ev else "",
        datetime.fromisoformat(ev[0]["created_at"].replace("Z", "+00:00")) if ev else None, "GitHub API repos/.../events (CreateEvent)")
except Exception as exc:  # events are kept for a limited time only
    add("independent_clock", "github_first_push_branch_created", f"not available ({type(exc).__name__})", None, "GitHub API events")

# ---------------------------------------------------------------------------------------------- Slurm record against the preregistration commit
jobs = []
with open(REPO / "final_checks" / "sacct_research_jobs.csv", encoding="utf-8", newline="") as fh:
    for r in csv.DictReader(fh):
        r["_sub"] = edt_naive(r["submit"]); r["_start"] = edt_naive(r["start"]); r["_end"] = edt_naive(r["end"])
        jobs.append(r)
study = [j for j in jobs if re.match(r"(lit|seaad|d1|d2|d3)", j["jobname"])]
before = [j for j in study if j["_sub"] and j["_sub"] < t_prereg]
add("before_preregistration_commit", "study_cluster_jobs_submitted_before_commit", str(len(before)), None,
    "final_checks/sacct_research_jobs.csv (sacct -X -P on ORCD; research jobs lit-, seaad-, d1-, d2-, d3; an unrelated animation-repo job `seer-install` is not in the table)")


def first(prefix, key="_sub"):
    c = [j for j in study if j["jobname"].startswith(prefix) and j[key]]
    return min(c, key=lambda j: j[key]) if c else None


for label, prefix, key in [("first_study_job_submitted", "", "_sub"),
                           ("first_bulk_data_download_job_started", "lit-dl", "_start"),
                           ("first_seaad_data_download_job_submitted", "seaad-dl", "_sub"),
                           ("first_seaad_data_inspection_job_started", "seaad-inspect", "_start"),
                           ("first_seaad_cna_job_started", "seaad-cna", "_start")]:
    j = first(prefix, key)
    secs = int((j[key] - t_prereg).total_seconds())
    add("after_preregistration_commit", label, f"job {j['jobid']} {j['jobname']} ({secs} s after the commit)", j[key], "final_checks/sacct_research_jobs.csv")

# ---------------------------------------------------------------------------------------------- what was looked at before the commit (planning, from the protocol)
src_proto = "protocol/06-research-project.md (pilot numbers, data inventory D2, Appendix A); preregistration.md section 0.1"
add("looked_at_before_preregistration_commit", "pilot_literature_and_genetics_counts",
    "free API calls: PubMed counts of Alzheimer's papers by year, TREM2 counts by year, a handful of PubTator3 co-mention counts with a rough Spearman correlation against genetic score among 32 top-genetic genes, Open Targets 26.09 target counts; they informed the literature-study expectations, which are therefore confirmatory of the pilot, not blind", None, src_proto)
add("looked_at_before_preregistration_commit", "seaad_bucket_listing_readme_and_csv_header",
    "listing of the SEA-AD S3 bucket, its README, and the header of the final-nuclei metadata CSV (HTTP range request)", None, src_proto)
add("looked_at_before_preregistration_commit", "seaad_metadata_summary",
    "one streaming pass over the SEA-AD MTG final-nuclei metadata CSV (not saved to disk): per-donor cohort summary (84 donors; APOE genotype counts; ADNC, Braak, LATE, Lewy body, cognitive status; 11 Severely Affected Donors; nuclei per subclass) and the median and range of the CPS scores (from the 55 KB CPS table)", None, src_proto)
add("not_done_before_preregistration_commit", "seaad_expression_data_download",
    "the CELLxGENE subclass h5ad files and the pseudobulk objects were not downloaded; none of the preregistered analyses had been run on study data (CNA, brick scores, programs, MAGMA, H4, or the literature metrics on the PubMed and PubTator3 corpus); the bulk literature downloads had not started", None,
    "final_checks/sacct_research_jobs.csv (no study job before the commit); the rows above")

# ---------------------------------------------------------------------------------------------- data sources
with open(REPO / "results" / "data_provenance.csv", encoding="utf-8", newline="") as fh:
    for r in csv.DictReader(fh):
        s = r["retrieved_utc"]
        dt = None
        m = re.match(r"^(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2})(?::(\d{2}))?Z?$", s)
        if m:
            dt = datetime.fromisoformat(f"{m.group(1)}T{m.group(2)}:{m.group(3) or '00'}+00:00")
        item = f"{r['source']}: {r['item']}"
        val = f"{r['version_or_size']} | {r['url_or_query']}" + (f" | {r['sha256_or_note']}" if r["sha256_or_note"] else "")
        note = "" if dt else f" (retrieved {s})"
        add("data_source_literature", item, val + note, dt, "results/data_provenance.csv (retrieved_utc as recorded there)")


def job(jobid):
    return next(j for j in jobs if j["jobid"] == jobid)


src_d0 = "Day-0 data report (research/12-seaad-day0.md in the project folder, section 2); final_checks/sacct_research_jobs.csv for the retrieval window"
j1, j2 = job("25618162"), job("25618163")
add("data_source_seaad", "SEA-AD MTG subclass h5ads (CZ CELLxGENE)",
    f"collection 1ca90a2d-2943-483d-b678-b809bf464c30; 13 subclass files, 29.84 GB; ORCD job 25618162 ended {j1['end']} EDT", j1["_start"], src_d0)
add("data_source_seaad", "SEA-AD AWS Open Data, release 2026-06-22",
    f"s3://sea-ad-single-cell-profiling: MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei_metadata.2026-06-22.csv (1,153,185,413 B); Multiregion_2026 Global_and_Local_CPS.20260105.csv, pseudobulk objects 2026-06-22, pertpy_compositional_modeling outputs (pertpy_summary_CPS_Local.20260622.csv); ORCD job 25618163 ended {j2['end']} EDT", j2["_start"], src_d0)
add("data_source_genetics", "Bellenguez et al. 2022 AD GWAS, harmonised summary statistics",
    "GCST90027158 (35379992-GCST90027158-MONDO_0004975.h.tsv.gz, 1,140,755,970 B, MD5 e2c8be73b5aa7698c5e8878ae607fe85); https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/harmonised/; same job 25618163", j2["_start"], src_d0)
add("data_source_genetics", "MAGMA",
    "v1.10 (linux static); 1000 Genomes European reference panel g1000_eur (build 37); gene locations NCBI37.3 and NCBI38 (MAGMA site files hosted at vu.data.surfsara.nl); same job 25618163", j2["_start"], src_d0)
add("data_source_software", "analysis software",
    "ORCD env `sc`: python 3.12, scanpy 1.12.4, anndata 0.13.4, cna 0.2.6, pertpy 1.4.0, pydeseq2 0.5.4, statsmodels 0.15.0, scikit-learn 1.9.1; literature env: see the Software row of results/data_provenance.csv", None, src_d0)

with open(OUT, "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["section", "item", "value", "time_edt", "time_utc", "source"])
    w.writerows(rows)
print(f"wrote {OUT} ({len(rows)} rows)")
