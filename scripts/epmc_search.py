#!/usr/bin/env python3
"""Europe PMC search helper (public REST API, no personal data sent)."""
import json, sys, urllib.parse, urllib.request

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

def search(q, n=8):
    url = BASE + "?" + urllib.parse.urlencode(
        {"query": q, "format": "json", "resultType": "core", "pageSize": n})
    req = urllib.request.Request(url, headers={"User-Agent": "literature-search-client"})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return json.load(fh)

if __name__ == "__main__":
    q = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    d = search(q, n)
    print(f"### {q}   -> {d.get('hitCount')} hits")
    for r in d.get("resultList", {}).get("result", []):
        ji = r.get("journalInfo", {}) or {}
        print(f"- {r.get('pubYear')} {(ji.get('journal') or {}).get('title','?')[:40]} | "
              f"PMID {r.get('pmid')} PMC {r.get('pmcid')} OA={r.get('isOpenAccess')}")
        print(f"    {(r.get('title') or '')[:150]}")
