#!/usr/bin/env bash
set -euo pipefail

MODE=${1:-p2p}
SCENARIO=${2:-biovault/tests/scenarios/syqure-distributed.yaml}

case "$MODE" in
  p2p|quic)
    export BV_SYFTBOX_HOTLINK=1
    export BV_SYFTBOX_HOTLINK_QUIC=1
    export BV_SYFTBOX_HOTLINK_QUIC_ONLY=1
    export BV_SYQURE_TRANSPORT=hotlink
    ;;
  ws-only|hotlink-ws)
    export BV_SYFTBOX_HOTLINK=1
    export BV_SYFTBOX_HOTLINK_QUIC=0
    export BV_SYFTBOX_HOTLINK_QUIC_ONLY=0
    export BV_SYQURE_TRANSPORT=hotlink
    ;;
  file)
    export BV_SYFTBOX_HOTLINK=0
    export BV_SYFTBOX_HOTLINK_QUIC=0
    export BV_SYFTBOX_HOTLINK_QUIC_ONLY=0
    export BV_SYQURE_TRANSPORT=file
    ;;
  *)
    echo "Usage: $0 {p2p|ws-only|file} [scenario]" >&2
    echo "Aliases: quic -> p2p, hotlink-ws -> ws-only" >&2
    exit 2
    ;;
 esac

echo "Mode: $MODE"
echo "Scenario: $SCENARIO"

# Use /usr/bin/time for consistent output formatting.
/usr/bin/time -p ./biovault/test-scenario.sh "$SCENARIO"
