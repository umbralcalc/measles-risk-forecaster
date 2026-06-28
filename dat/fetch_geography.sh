#!/usr/bin/env bash
#
# Fetch England Upper Tier Local Authority boundaries from the ONS Open
# Geography Portal (gating check #4 — the spatial layer).
#
# Layer: Counties and Unitary Authorities (December 2024) Boundaries UK BGC.
# This is the geography that matches the current COVER UTLA set: metropolitan
# districts are individual (153 England areas), unlike the older "Upper Tier
# Local Authorities" layer which groups them into metropolitan counties (122).
# The 153-vs-152 difference against COVER is the documented small-LA combination
# (e.g. City of London, Isles of Scilly) — reconciled at the join, not here.
#
# Output: dat/utla_boundaries.geojson (CTYUA24CD, CTYUA24NM + geometry, WGS84).
#
# Licence: Office for National Statistics, Open Government Licence v3.0.
# Contains OS data © Crown copyright and database right 2024. (See SOURCES.md §4.)
#
# Usage: ./dat/fetch_geography.sh

set -eo pipefail

DATA_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT="${DATA_DIR}/utla_boundaries.geojson"

BASE="https://services1.arcgis.com/ESMARspQHYMw9BZ9/arcgis/rest/services"
LAYER="Counties_and_Unitary_Authorities_December_2024_Boundaries_UK_BGC"

echo "Fetching England UTLA (CUA Dec 2024) boundaries from ONS Open Geography Portal..."
curl -sf --max-time 120 -G "${BASE}/${LAYER}/FeatureServer/0/query" \
  --data-urlencode "where=CTYUA24CD LIKE 'E%'" \
  --data-urlencode "outFields=CTYUA24CD,CTYUA24NM" \
  --data-urlencode "returnGeometry=true" \
  --data-urlencode "outSR=4326" \
  --data-urlencode "f=geojson" \
  -o "${OUT}"

# Sanity: report feature count and flag if the transfer limit was hit.
python3 - "${OUT}" <<'PY'
import json, sys
with open(sys.argv[1]) as f:
    gj = json.load(f)
feats = gj.get("features", [])
print(f"  -> {sys.argv[1].split('/')[-1]}: {len(feats)} features")
if gj.get("exceededTransferLimit"):
    print("  WARNING: exceededTransferLimit=true — pagination needed, data is incomplete!", file=sys.stderr)
    sys.exit(1)
if len(feats) != 153:
    print(f"  WARNING: expected 153 England areas, got {len(feats)}", file=sys.stderr)
PY

echo "Done."
