#!/bin/bash
set -euo pipefail

# Generate rsids and counts based on the current datasite
# BV_CURRENT_DATASITE is set by the flow runner

case "${BV_CURRENT_DATASITE:-}" in
  client1@sandbox.local)
    printf "rs1\nrs2\n" > rsids.txt
    echo '{"rs1": 3, "rs2": 1}' > counts.json
    ;;
  client2@sandbox.local)
    printf "rs2\nrs3\nrs4\n" > rsids.txt
    echo '{"rs2": 2, "rs3": 4, "rs4": 5}' > counts.json
    ;;
  *)
    echo "Unknown datasite: ${BV_CURRENT_DATASITE:-unset}" >&2
    exit 1
    ;;
esac

echo "Generated rsids and counts for ${BV_CURRENT_DATASITE}"
