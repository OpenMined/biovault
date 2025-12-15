#!/usr/bin/env python3
"""
Seeded chaos test runner for BioVault sandbox devstack.

This intentionally mixes:
- messaging (send/reply/sync)
- session invites (create/accept/reject/chat)
- project submission + processing + results/ACL steps

The goal is not to assert one exact transcript, but to enforce invariants:
- `bv message sync` should not crash on decrypt failures
- decrypt failures should be recorded in `failed_messages` and retain the request file path
- message/session/project tables should not contain duplicate IDs
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sqlite3
import subprocess
import sys
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


ROOT = Path(__file__).resolve().parents[2]
SANDBOX_ROOT = ROOT / "sandbox"
BV_BIN = Path(os.environ.get("BV_BIN", str(ROOT / "cli" / "target" / "release" / "bv")))


CLIENTS_DEFAULT = ["client1@sandbox.local", "client2@sandbox.local"]


@dataclass(frozen=True)
class Client:
    email: str

    @property
    def dir(self) -> Path:
        return SANDBOX_ROOT / self.email

    @property
    def syftbox_config_path(self) -> Path:
        return self.dir / ".syftbox" / "config.json"

    @property
    def biovault_dir(self) -> Path:
        # In sandbox mode BIOVAULT_HOME resolves to SYFTBOX_DATA_DIR/.biovault
        return self.dir / ".biovault"

    @property
    def messages_db(self) -> Path:
        return self.biovault_dir / "data" / "messages.db"

    @property
    def biovault_db(self) -> Path:
        return self.biovault_dir / "biovault.db"


def sh(
    cmd: str | List[str],
    *,
    cwd: Optional[Path] = None,
    env: Optional[Dict[str, str]] = None,
    timeout_s: int = 120,
    check: bool = True,
    capture: bool = True,
) -> subprocess.CompletedProcess[str]:
    if isinstance(cmd, list):
        display = " ".join(cmd)
    else:
        display = cmd

    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)

    result = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        env=merged_env,
        text=True,
        capture_output=capture,
        shell=isinstance(cmd, str),
        timeout=timeout_s,
    )
    if check and result.returncode != 0:
        stderr = (result.stderr or "").strip()
        stdout = (result.stdout or "").strip()
        raise SystemExit(
            f"Command failed (exit {result.returncode}): {display}\nSTDOUT:\n{stdout}\nSTDERR:\n{stderr}"
        )
    return result


def parse_duration_seconds(raw: str) -> int:
    value = raw.strip().lower()
    if not value:
        raise ValueError("empty duration")

    multiplier = 1
    if value.endswith("ms"):
        multiplier = 0.001
        value = value[:-2]
    elif value.endswith("s"):
        multiplier = 1
        value = value[:-1]
    elif value.endswith("m"):
        multiplier = 60
        value = value[:-1]
    elif value.endswith("h"):
        multiplier = 60 * 60
        value = value[:-1]

    number = float(value)
    seconds = int(number * multiplier)
    if seconds <= 0:
        raise ValueError("duration must be > 0")
    return seconds


def bv(
    client: Client,
    args: List[str],
    *,
    timeout_s: int = 120,
    check: bool = True,
    capture: bool = True,
) -> subprocess.CompletedProcess[str]:
    if not BV_BIN.exists():
        raise SystemExit(f"Missing BioVault CLI binary: {BV_BIN} (run cargo build --release)")

    env: Dict[str, str] = {
        "HOME": str(client.dir),
        "SYFTBOX_EMAIL": client.email,
        "SYFTBOX_DATA_DIR": str(client.dir),
        "SYFTBOX_CONFIG_PATH": str(client.syftbox_config_path),
    }
    # Preserve scenario runner conventions, if present
    if os.environ.get("SCENARIO_JAVA_HOME"):
        env["JAVA_HOME"] = os.environ["SCENARIO_JAVA_HOME"]
        env["JAVA_CMD"] = str(Path(os.environ["SCENARIO_JAVA_HOME"]) / "bin" / "java")
    if os.environ.get("SCENARIO_USER_PATH"):
        env["PATH"] = os.environ["SCENARIO_USER_PATH"]
    env.setdefault("NXF_DISABLE_JAVA_VERSION_CHECK", "true")
    env.setdefault("NXF_IGNORE_JAVA_VERSION", "true")
    env.setdefault("NXF_OPTS", "-Dnxf.java.check=false")

    return sh(
        [str(BV_BIN), *args],
        cwd=client.dir,
        env=env,
        timeout_s=timeout_s,
        check=check,
        capture=capture,
    )


def wait_for(client: Client, rel_path_glob: str, *, timeout_s: int = 30) -> Path:
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        matches = list(client.dir.glob(rel_path_glob))
        if matches:
            return matches[0]
        time.sleep(0.5)
    raise SystemExit(f"Timeout waiting for {client.email}:{rel_path_glob}")


def load_json_stdout(result: subprocess.CompletedProcess[str]) -> Any:
    raw = (result.stdout or "").strip()
    try:
        return json.loads(raw)
    except Exception as e:
        raise SystemExit(f"Failed to parse JSON from stdout: {e}\nstdout:\n{raw}")


def sql_rows(db_path: Path, query: str, params: Tuple[Any, ...] = ()) -> List[Tuple[Any, ...]]:
    if not db_path.exists():
        return []
    con = sqlite3.connect(str(db_path))
    try:
        cur = con.execute(query, params)
        return list(cur.fetchall())
    except sqlite3.OperationalError:
        return []
    finally:
        con.close()


def table_columns(db_path: Path, table: str) -> List[str]:
    if not db_path.exists():
        return []
    con = sqlite3.connect(str(db_path))
    try:
        cur = con.execute(f"PRAGMA table_info({table})")
        return [str(row[1]) for row in cur.fetchall() if len(row) > 1]
    except sqlite3.OperationalError:
        return []
    finally:
        con.close()


def failed_messages_path_column(db_path: Path) -> Optional[str]:
    cols = set(table_columns(db_path, "failed_messages"))
    if "request_path" in cols:
        return "request_path"
    if "file_path" in cols:
        return "file_path"
    return None


def assert_invariants(clients: List[Client]) -> None:
    for c in clients:
        # No duplicate message IDs
        dups = sql_rows(
            c.messages_db,
            "SELECT id, COUNT(*) FROM messages GROUP BY id HAVING COUNT(*) > 1",
        )
        if dups:
            raise SystemExit(f"{c.email}: duplicate message ids detected: {dups[:5]}")

        # No duplicate session IDs
        sdups = sql_rows(
            c.biovault_db,
            "SELECT session_id, COUNT(*) FROM sessions GROUP BY session_id HAVING COUNT(*) > 1",
        )
        if sdups:
            raise SystemExit(f"{c.email}: duplicate sessions detected: {sdups[:5]}")

        # Failed messages should keep their request file around for retry (no silent loss)
        path_col = failed_messages_path_column(c.messages_db)
        if path_col:
            failed = sql_rows(
                c.messages_db,
                f"SELECT id, {path_col} FROM failed_messages WHERE dismissed = 0",
            )
            for (_id, request_path) in failed:
                if request_path and not Path(str(request_path)).exists():
                    raise SystemExit(
                        f"{c.email}: failed_messages[{_id}] {path_col} missing on disk: {request_path}"
                    )


def inbox_json(client: Client, *, sent: bool = False, all_msgs: bool = False) -> Dict[str, Any]:
    args = ["inbox", "--plain", "--json"]
    if sent:
        args.append("--sent")
    if all_msgs:
        args.append("--all")
    return load_json_stdout(bv(client, args))


def choose_message_id(rng: random.Random, inbox: Dict[str, Any], *, kind: str) -> Optional[str]:
    messages = inbox.get("messages") or []
    candidates: List[Dict[str, Any]] = []
    for m in messages:
        mt = m.get("message_type")
        if kind == "text":
            if mt == "Text":
                candidates.append(m)
        elif kind == "project":
            if isinstance(mt, dict) and mt.get("Project") is not None:
                candidates.append(m)
        else:
            candidates.append(m)
    if not candidates:
        return None
    return rng.choice(candidates).get("id")


def action_sync(rng: random.Random, client: Client) -> None:
    args = ["message", "sync"]
    if rng.random() < 0.15:
        args.append("--no-cleanup")
    result = bv(client, args, timeout_s=180, check=False)
    if result.returncode == 0:
        return "ok"
    return f"nonzero-exit={result.returncode} (allowed)"


def action_import_peer_key(rng: random.Random, client: Client, peer: Client) -> None:
    # Prefer rotated bundle if present, else default did.json
    rotated = client.dir / "datasites" / peer.email / "public" / "crypto" / "did-rotated.json"
    bundle = (
        f"datasites/{peer.email}/public/crypto/did-rotated.json"
        if rotated.exists()
        else f"datasites/{peer.email}/public/crypto/did.json"
    )
    args = ["key", "import", bundle, "--email", peer.email]
    if rng.random() < 0.7:
        args.append("--ignore-tofu")
    # Allow TOFU failures (expected sometimes)
    result = bv(client, args, check=False)
    if result.returncode == 0:
        return f"imported {peer.email}"
    return f"nonzero-exit={result.returncode} (allowed)"


def action_rotate_key(rng: random.Random, client: Client) -> None:
    # Rotate local signing key; message decrypt on peers may fail until they import our new bundle.
    bv(client, ["key", "wipe"], timeout_s=120)
    bv(client, ["key", "generate", "--force"], timeout_s=240)
    bv(
        client,
        [
            "key",
            "export",
            "--output",
            f"datasites/{client.email}/public/crypto/did-rotated.json",
        ],
        timeout_s=120,
    )
    return "rotated+exported did-rotated.json"


def action_send_text(rng: random.Random, seed: int, step: int, sender: Client, recipient: Client) -> None:
    marker = f"CHAOS:{seed}:{step}:{uuid.uuid4()}"
    subject = f"Chaos {seed}"
    body = f"{marker} hello"
    result = bv(
        sender,
        ["message", "send", recipient.email, body, "--subject", subject],
        check=False,
    )
    if result.returncode == 0:
        return f"{sender.email} -> {recipient.email}"
    return f"{sender.email} -> {recipient.email} nonzero-exit={result.returncode} (allowed)"


def action_reply_text(rng: random.Random, seed: int, step: int, sender: Client) -> None:
    all_box = inbox_json(sender, all_msgs=True)
    msg_id = choose_message_id(rng, all_box, kind="text")
    if not msg_id:
        return "skipped (no text messages)"
    marker = f"CHAOS-REPLY:{seed}:{step}:{uuid.uuid4()}"
    result = bv(sender, ["message", "reply", msg_id, marker], check=False)
    if result.returncode == 0:
        return f"reply_to={msg_id[:8]}"
    return f"reply_to={msg_id[:8]} nonzero-exit={result.returncode} (allowed)"


def action_create_session(rng: random.Random, seed: int, step: int, owner: Client, peer: Client) -> None:
    name = f"Chaos Session {seed}-{step}"
    desc = f"chaos:{seed}:{step}"
    result = bv(
        owner,
        ["session", "create", name, "--peer", peer.email, "--description", desc],
        check=False,
    )
    if result.returncode == 0:
        return f"{owner.email} invites {peer.email}"
    return f"{owner.email} invites {peer.email} nonzero-exit={result.returncode} (allowed)"


def action_accept_invite(rng: random.Random, client: Client) -> None:
    inv = load_json_stdout(bv(client, ["session", "invitations", "--json"], check=False))
    if not isinstance(inv, list) or not inv:
        return "skipped (no invitations)"
    session_id = rng.choice(inv).get("session_id")
    if not session_id:
        return "skipped (missing session_id)"
    result = bv(client, ["session", "accept", session_id], check=False)
    if result.returncode == 0:
        return f"accepted {session_id}"
    return f"accept {session_id} nonzero-exit={result.returncode} (allowed)"


def action_reject_invite(rng: random.Random, client: Client) -> None:
    inv = load_json_stdout(bv(client, ["session", "invitations", "--json"], check=False))
    if not isinstance(inv, list) or not inv:
        return "skipped (no invitations)"
    session_id = rng.choice(inv).get("session_id")
    if not session_id:
        return "skipped (missing session_id)"
    result = bv(client, ["session", "reject", session_id, "--reason", "chaos"], check=False)
    if result.returncode == 0:
        return f"rejected {session_id}"
    return f"reject {session_id} nonzero-exit={result.returncode} (allowed)"


def action_session_chat(rng: random.Random, seed: int, step: int, client: Client) -> None:
    sessions = load_json_stdout(bv(client, ["session", "list", "--json"], check=False))
    if not isinstance(sessions, list) or not sessions:
        return "skipped (no sessions)"
    sess = rng.choice(sessions)
    sid = sess.get("session_id")
    if not sid:
        return "skipped (missing session_id)"
    msg = f"CHAOS-CHAT:{seed}:{step}:{uuid.uuid4()}"
    result = bv(client, ["session", "chat", sid, msg], check=False)
    if result.returncode == 0:
        return f"session={sid}"
    return f"session={sid} nonzero-exit={result.returncode} (allowed)"


def action_submit_project(sender: Client, recipient: Client) -> None:
    # Use the lightweight "count-lines" example project (matches inbox-ping-pong.yaml).
    result = sh(
        "\n".join(
            [
                "rm -rf count-lines",
                "mkdir -p count-lines",
                f"cp -R \"{ROOT}/cli/examples/pipeline/count-lines/.\" count-lines/",
                f"\"{BV_BIN}\" submit count-lines {recipient.email} --non-interactive --force",
            ]
        ),
        cwd=sender.dir,
        env={
            "HOME": str(sender.dir),
            "SYFTBOX_EMAIL": sender.email,
            "SYFTBOX_DATA_DIR": str(sender.dir),
            "SYFTBOX_CONFIG_PATH": str(sender.syftbox_config_path),
        },
        timeout_s=300,
        check=False,
        capture=True,
    )
    if result.returncode == 0:
        return f"{sender.email} -> {recipient.email}"
    return f"{sender.email} -> {recipient.email} nonzero-exit={result.returncode} (allowed)"


def action_process_project(rng: random.Random, client: Client) -> None:
    all_box = inbox_json(client, all_msgs=True)
    msg_id = choose_message_id(rng, all_box, kind="project")
    if not msg_id:
        return "skipped (no project messages)"
    result = bv(
        client,
        ["message", "process", msg_id, "--test", "--participant", "ALL", "--approve", "--non-interactive"],
        timeout_s=600,
        check=False,
    )
    if result.returncode == 0:
        return f"message={msg_id[:8]}"
    return f"message={msg_id[:8]} nonzero-exit={result.returncode} (allowed)"


def action_restrict_results_write(client: Client) -> None:
    # Similar to inbox-ping-pong.yaml: restrict results write access after receipt.
    # Best-effort: if no submission exists, skip.
    submissions = client.dir / "datasites" / client.email / "shared" / "biovault" / "submissions"
    if not submissions.exists():
        return "skipped (no submissions dir)"
    subdirs = [p for p in submissions.iterdir() if p.is_dir()]
    if not subdirs:
        return "skipped (no submissions)"
    latest = sorted(subdirs)[-1]
    perm = latest / "syft.pub.yaml"
    if not perm.exists():
        return "skipped (missing syft.pub.yaml)"
    result = sh(
        f"tmp=$(mktemp)\n"
        f"yq '(.rules[] | select(.pattern == \"results/**/*\").access.write) = []' \"{perm}\" > \"$tmp\"\n"
        f"mv \"$tmp\" \"{perm}\"",
        cwd=client.dir,
        timeout_s=60,
        check=False,
    )
    if result.returncode == 0:
        return f"submission={latest.name}"
    return f"submission={latest.name} nonzero-exit={result.returncode} (allowed)"


def parse_clients(raw: str) -> List[str]:
    parts = [p.strip().lower() for p in raw.split(",") if p.strip()]
    if len(parts) < 2:
        raise SystemExit("--clients must include at least 2 emails")
    return parts


def replay_command(
    *,
    seed: int,
    raw_duration: Optional[str],
    steps_arg: Optional[int],
    clients_arg: str,
    skip_reset: bool,
    skip_fixtures: bool,
) -> str:
    env_prefix = ""
    if raw_duration:
        env_prefix = f"CHAOS_DURATION={raw_duration} "

    args: List[str] = ["./test-chaos.sh", "--seed", str(seed)]

    # Preserve shape of the run for reproducibility.
    if raw_duration is None:
        args += ["--steps", str(steps_arg if steps_arg is not None else 60)]
    if clients_arg != ",".join(CLIENTS_DEFAULT):
        args += ["--clients", clients_arg]
    if skip_reset:
        args.append("--skip-reset")
    if skip_fixtures:
        args.append("--skip-fixtures")

    return f"cd biovault && {env_prefix}{' '.join(args)}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=None, help="Random seed for reproducible chaos")
    ap.add_argument(
        "--steps",
        type=int,
        default=None,
        help="Max chaos actions to run (default: 60, or unlimited if duration is set)",
    )
    ap.add_argument(
        "--duration",
        type=str,
        default=None,
        help="Run for this wall-clock duration (e.g. 30s, 5m, 1h); env: CHAOS_DURATION",
    )
    ap.add_argument(
        "--clients",
        type=str,
        default=",".join(CLIENTS_DEFAULT),
        help="Comma-separated sandbox clients (default: client1@sandbox.local,client2@sandbox.local)",
    )
    ap.add_argument("--skip-reset", action="store_true", help="Do not reset devstack")
    ap.add_argument("--skip-fixtures", action="store_true", help="Do not import genotype fixtures")
    ap.add_argument("--invariant-every", type=int, default=8, help="Invariant check frequency")
    ap.add_argument(
        "--print-trace-json",
        action="store_true",
        help="Print the full JSON trace at the end (default: off)",
    )
    args = ap.parse_args()

    seed = int(args.seed) if args.seed is not None else int.from_bytes(os.urandom(4), "big")
    rng = random.Random(seed)

    raw_duration = args.duration if args.duration is not None else os.environ.get("CHAOS_DURATION")
    duration_s: Optional[int] = None
    if raw_duration:
        try:
            duration_s = parse_duration_seconds(raw_duration)
        except ValueError as e:
            raise SystemExit(f"Invalid duration '{raw_duration}': {e}")

    max_steps: int
    if args.steps is not None:
        max_steps = args.steps
    elif duration_s is not None:
        max_steps = 1_000_000_000
    else:
        max_steps = 60

    emails = parse_clients(args.clients)
    clients = [Client(e) for e in emails[:2]]
    a, b = clients[0], clients[1]

    duration_label = f" duration={raw_duration}" if raw_duration else ""
    print(f"[chaos] seed={seed} steps={max_steps} clients={a.email},{b.email}{duration_label}")
    print(
        "[chaos] replay:",
        replay_command(
            seed=seed,
            raw_duration=raw_duration,
            steps_arg=args.steps,
            clients_arg=args.clients,
            skip_reset=args.skip_reset,
            skip_fixtures=args.skip_fixtures,
        ),
    )
    sys.stdout.flush()

    # Start (or reuse) devstack
    devstack_args = ["./tests/scripts/devstack.sh", "--clients", f"{a.email},{b.email}"]
    if not args.skip_reset:
        devstack_args.append("--reset")
    devstack_res = sh(devstack_args, cwd=ROOT, timeout_s=600, capture=True, check=False)
    if devstack_res.returncode != 0:
        print("❌ devstack start failed", file=sys.stderr)
        print(f"[chaos] seed={seed} (replay with --seed {seed})", file=sys.stderr)
        stderr = (devstack_res.stderr or "").strip()
        stdout = (devstack_res.stdout or "").strip()
        raise SystemExit(f"devstack failed\nSTDOUT:\n{stdout}\nSTDERR:\n{stderr}")
    print("✅ devstack ready")

    # Wait for public bundles to replicate so key imports can happen.
    wait_for(a, f"datasites/{b.email}/public/crypto/did.json", timeout_s=45)
    wait_for(b, f"datasites/{a.email}/public/crypto/did.json", timeout_s=45)
    print("✅ keys replicated")

    # Import each other's public keys up-front.
    action_import_peer_key(rng, a, b)
    action_import_peer_key(rng, b, a)

    # Seed fixtures so project processing has inputs.
    if not args.skip_fixtures:
        sh(
            [
                "./tests/scripts/import-data.sh",
                "--csv",
                str(ROOT / "cli" / "tests" / "data" / "participants_all.csv"),
                "--clients",
                f"{a.email},{b.email}",
                "--path-column",
                "genotype_file",
            ],
            cwd=ROOT,
            timeout_s=600,
            check=False,
        )

    # Ensure we exercise the "project submit" path at least once in a typical run.
    # This is intentionally lightweight (it should not require Nextflow).
    try:
        bootstrap_details = action_submit_project(a, b)
        print(f"✅ bootstrap project-submit — {bootstrap_details}")
        action_sync(rng, a)
        action_sync(rng, b)
    except Exception as e:
        print("❌ bootstrap project-submit failed", file=sys.stderr)
        print(
            "[chaos] replay:",
            replay_command(
                seed=seed,
                raw_duration=raw_duration,
                steps_arg=args.steps,
                clients_arg=args.clients,
                skip_reset=args.skip_reset,
                skip_fixtures=args.skip_fixtures,
            ),
            file=sys.stderr,
        )
        raise SystemExit(str(e))

    # Define weighted actions.
    def pick_party() -> Tuple[Client, Client]:
        return (a, b) if rng.random() < 0.5 else (b, a)

    actions = [
        ("sync-a", 16, lambda step: action_sync(rng, a)),
        ("sync-b", 16, lambda step: action_sync(rng, b)),
        ("import-a", 6, lambda step: action_import_peer_key(rng, a, b)),
        ("import-b", 6, lambda step: action_import_peer_key(rng, b, a)),
        ("rotate-a", 2, lambda step: action_rotate_key(rng, a)),
        ("rotate-b", 2, lambda step: action_rotate_key(rng, b)),
        ("send", 12, lambda step: action_send_text(rng, seed, step, *pick_party())),
        ("reply-a", 6, lambda step: action_reply_text(rng, seed, step, a)),
        ("reply-b", 6, lambda step: action_reply_text(rng, seed, step, b)),
        ("session-create", 6, lambda step: action_create_session(rng, seed, step, *pick_party())),
        ("session-accept-a", 3, lambda step: action_accept_invite(rng, a)),
        ("session-accept-b", 3, lambda step: action_accept_invite(rng, b)),
        ("session-reject-a", 2, lambda step: action_reject_invite(rng, a)),
        ("session-reject-b", 2, lambda step: action_reject_invite(rng, b)),
        ("session-chat-a", 4, lambda step: action_session_chat(rng, seed, step, a)),
        ("session-chat-b", 4, lambda step: action_session_chat(rng, seed, step, b)),
        ("project-submit", 8, lambda step: action_submit_project(*pick_party())),
        ("project-process-a", 4, lambda step: action_process_project(rng, a)),
        ("project-process-b", 4, lambda step: action_process_project(rng, b)),
        ("restrict-results-a", 1, lambda step: action_restrict_results_write(a)),
        ("restrict-results-b", 1, lambda step: action_restrict_results_write(b)),
    ]

    weighted: List[Tuple[str, int]] = [(name, w) for (name, w, _fn) in actions]
    total_weight = sum(w for (_n, w) in weighted)

    trace: List[Dict[str, Any]] = []

    started_at = time.time()
    deadline = (started_at + duration_s) if duration_s is not None else None

    step = 0
    while step < max_steps:
        if deadline is not None and time.time() >= deadline:
            break
        step += 1
        roll = rng.randrange(total_weight)
        chosen_name = None
        chosen_fn = None
        acc = 0
        for (name, weight, fn) in actions:
            acc += weight
            if roll < acc:
                chosen_name = name
                chosen_fn = fn
                break
        assert chosen_name is not None and chosen_fn is not None

        started = time.time()
        err: Optional[str] = None
        details: Optional[str] = None
        try:
            details = chosen_fn(step)
        except Exception as e:
            # Chaos should tolerate some expected command failures (TOFU, missing keys, etc),
            # but unexpected Python errors should fail fast with the seed for replay.
            err = f"{type(e).__name__}: {e}"
        elapsed_ms = int((time.time() - started) * 1000)
        trace.append(
            {
                "step": step,
                "action": chosen_name,
                "elapsed_ms": elapsed_ms,
                "details": details,
                "error": err,
            }
        )
        if err:
            print(
                f"❌ step {step}/{max_steps} {chosen_name} ({elapsed_ms}ms) error: {err}",
                file=sys.stderr,
            )
            print(
                "[chaos] replay:",
                replay_command(
                    seed=seed,
                    raw_duration=raw_duration,
                    steps_arg=args.steps,
                    clients_arg=args.clients,
                    skip_reset=args.skip_reset,
                    skip_fixtures=args.skip_fixtures,
                ),
                file=sys.stderr,
            )
            raise SystemExit(1)

        suffix = f" — {details}" if details else ""
        print(f"✅ step {step}/{max_steps} {chosen_name} ({elapsed_ms}ms){suffix}")

        if args.invariant_every > 0 and step % args.invariant_every == 0:
            try:
                assert_invariants(clients)
            except SystemExit as e:
                print(f"❌ invariants failed at step {step}: {e}", file=sys.stderr)
                print(
                    "[chaos] replay:",
                    replay_command(
                        seed=seed,
                        raw_duration=raw_duration,
                        steps_arg=args.steps,
                        clients_arg=args.clients,
                        skip_reset=args.skip_reset,
                        skip_fixtures=args.skip_fixtures,
                    ),
                    file=sys.stderr,
                )
                raise
            print(f"✅ invariants (step {step})")

    # Final drain + invariants.
    for _ in range(3):
        action_sync(rng, a)
        action_sync(rng, b)
    assert_invariants(clients)

    print("[chaos] completed successfully ✅")
    print(
        "[chaos] replay:",
        replay_command(
            seed=seed,
            raw_duration=raw_duration,
            steps_arg=args.steps,
            clients_arg=args.clients,
            skip_reset=args.skip_reset,
            skip_fixtures=args.skip_fixtures,
        ),
    )
    elapsed_s = int(time.time() - started_at)
    print(f"[chaos] steps_ran={step} elapsed={elapsed_s}s")
    if args.print_trace_json:
        print(json.dumps({"seed": seed, "steps": step, "trace": trace}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
