#!/usr/bin/env python3
"""Simple YAML-driven scenario runner for sandbox datasites."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Optional

ROOT = Path(__file__).resolve().parent.parent
SANDBOX_ROOT = ROOT / "sandbox"
JAVA_HOME_OVERRIDE: Optional[str] = None


def detect_java_home() -> Optional[str]:
    """Return a Java 17-24 home if available, preferring SCENARIO_JAVA_HOME."""
    env_java = os.environ.get("SCENARIO_JAVA_HOME") or os.environ.get("JAVA_HOME")
    if env_java and is_supported_java(env_java):
        return env_java

    helper = Path("/usr/libexec/java_home")
    if helper.exists():
        try:
            result = subprocess.run(
                [str(helper), "-v", "17-24"],
                check=True,
                capture_output=True,
                text=True,
            )
            candidate = result.stdout.strip()
            if candidate:
                return candidate
        except Exception:
            pass
    return None


def is_supported_java(java_home: str) -> bool:
    """Check whether the given JAVA_HOME points to a Java 17-24 runtime."""
    java_bin = Path(java_home) / "bin" / "java"
    if not java_bin.exists():
        return False
    try:
        result = subprocess.run(
            [str(java_bin), "-version"],
            capture_output=True,
            text=True,
        )
        output = result.stderr or result.stdout
        # Parse version like 'java version "17.0.12"' or 'openjdk version "25.0.1"'
        for line in output.splitlines():
            if "version" in line and '"' in line:
                try:
                    version_str = line.split('"')[1]
                    major = int(version_str.split(".")[0])
                    return 17 <= major <= 24
                except Exception:
                    continue
    except Exception:
        return False
    return False


JAVA_HOME_OVERRIDE = detect_java_home()


def load_scenario(path: Path) -> Dict[str, Any]:
    text = path.read_text()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        try:
            import yaml  # type: ignore
        except ImportError as exc:  # pragma: no cover - fallback message
            return minimal_yaml_load(text, path)
        return yaml.safe_load(text)  # type: ignore[arg-type]


def minimal_yaml_load(text: str, path: Path) -> Dict[str, Any]:
    """
    Minimal YAML loader for our scenario files.

    Supports a small, indentation-based subset of YAML used in this repo:
    - mappings (key: value, key: <nested>)
    - lists (- item, - key: value, - <nested mapping>)
    - scalars: strings (quoted/unquoted), ints, bools, null
    - block scalars: | and > (kept as raw text; > folds newlines to spaces)

    This intentionally rejects advanced YAML features (anchors, tags, complex
    multiline scalars, flow-style collections) to avoid silently misparsing.
    """

    class ParseError(Exception):
        pass

    def is_blank_or_comment(line: str) -> bool:
        stripped = line.strip()
        return not stripped or stripped.startswith("#")

    def indent_of(line: str) -> int:
        return len(line) - len(line.lstrip(" "))

    def strip_inline_comment(value: str) -> str:
        # For our subset, treat " # ..." as a comment delimiter when not quoted.
        # If a user needs a literal '#', they can quote the value.
        v = value.rstrip()
        if not v:
            return v
        if v[0] in ("'", '"'):
            return v
        hash_idx = v.find(" #")
        if hash_idx != -1:
            return v[:hash_idx].rstrip()
        return v

    def parse_scalar(raw: str) -> Any:
        v = strip_inline_comment(raw.strip())
        if v in ("null", "NULL", "~", ""):
            return None
        if v in ("true", "True", "TRUE"):
            return True
        if v in ("false", "False", "FALSE"):
            return False
        # int
        if v.isdigit() or (v.startswith("-") and v[1:].isdigit()):
            try:
                return int(v, 10)
            except Exception:
                pass
        # quoted string
        if len(v) >= 2 and ((v[0] == v[-1] == "'") or (v[0] == v[-1] == '"')):
            inner = v[1:-1]
            if v[0] == '"':
                inner = (
                    inner.replace(r"\\", "\\")
                    .replace(r"\"", '"')
                    .replace(r"\n", "\n")
                    .replace(r"\t", "\t")
                )
            return inner
        return v

    lines = text.splitlines()
    n = len(lines)

    def next_nonempty(idx: int) -> int:
        while idx < n and is_blank_or_comment(lines[idx]):
            idx += 1
        return idx

    def parse_block_scalar(start_idx: int, base_indent: int, style: str) -> tuple[str, int]:
        idx = start_idx
        buf: list[str] = []
        idx = next_nonempty(idx)
        while idx < n:
            line = lines[idx]
            if is_blank_or_comment(line):
                buf.append("")
                idx += 1
                continue
            ind = indent_of(line)
            if ind < base_indent:
                break
            buf.append(line[base_indent:])
            idx += 1
        if style == ">":
            # simple fold: join non-empty lines with spaces, preserve blank lines as newlines
            out_parts: list[str] = []
            paragraph: list[str] = []
            for l in buf:
                if l == "":
                    if paragraph:
                        out_parts.append(" ".join(paragraph))
                        paragraph = []
                    out_parts.append("")
                else:
                    paragraph.append(l)
            if paragraph:
                out_parts.append(" ".join(paragraph))
            return "\n".join(out_parts).rstrip("\n"), idx
        return "\n".join(buf).rstrip("\n"), idx

    def parse_node(idx: int, current_indent: int) -> tuple[Any, int]:
        idx = next_nonempty(idx)
        if idx >= n:
            return {}, idx

        # Decide whether we're parsing a list or a mapping at this level.
        line = lines[idx]
        if indent_of(line) < current_indent:
            return {}, idx

        if line.lstrip(" ").startswith("- "):
            items: list[Any] = []
            while True:
                idx = next_nonempty(idx)
                if idx >= n:
                    break
                line = lines[idx]
                ind = indent_of(line)
                if ind < current_indent:
                    break
                if not line.lstrip(" ").startswith("- "):
                    break
                if ind != current_indent:
                    raise ParseError(
                        f"{path}: list item has unexpected indent (expected {current_indent}, got {ind})"
                    )

                after_dash = line.lstrip(" ")[2:]
                if after_dash.strip() == "":
                    # "-\n  key: value" style
                    idx += 1
                    child, idx = parse_node(idx, current_indent + 2)
                    items.append(child)
                    continue

                if ":" in after_dash:
                    k, rest = after_dash.split(":", 1)
                    key = k.strip()
                    rest = rest.lstrip()
                    item_map: Dict[str, Any] = {}
                    if rest in ("|", ">"):
                        idx += 1
                        scalar, idx = parse_block_scalar(idx, current_indent + 4, rest)
                        item_map[key] = scalar
                        # Allow additional sibling keys for this list item (common in scenarios)
                        next_idx = next_nonempty(idx)
                        if (
                            next_idx < n
                            and indent_of(lines[next_idx]) == current_indent + 2
                            and not lines[next_idx].lstrip(" ").startswith("- ")
                        ):
                            extra, idx = parse_node(idx, current_indent + 2)
                            if not isinstance(extra, dict):
                                raise ParseError(f"{path}: expected mapping continuation at line {next_idx + 1}")
                            item_map.update(extra)
                        items.append(item_map)
                        continue
                    if rest == "":
                        idx += 1
                        child, idx = parse_node(idx, current_indent + 4)
                        item_map[key] = child
                        next_idx = next_nonempty(idx)
                        if (
                            next_idx < n
                            and indent_of(lines[next_idx]) == current_indent + 2
                            and not lines[next_idx].lstrip(" ").startswith("- ")
                        ):
                            extra, idx = parse_node(idx, current_indent + 2)
                            if not isinstance(extra, dict):
                                raise ParseError(f"{path}: expected mapping continuation at line {next_idx + 1}")
                            item_map.update(extra)
                        items.append(item_map)
                        continue
                    idx += 1
                    item_map[key] = parse_scalar(rest)
                    next_idx = next_nonempty(idx)
                    if (
                        next_idx < n
                        and indent_of(lines[next_idx]) == current_indent + 2
                        and not lines[next_idx].lstrip(" ").startswith("- ")
                    ):
                        extra, idx = parse_node(idx, current_indent + 2)
                        if not isinstance(extra, dict):
                            raise ParseError(f"{path}: expected mapping continuation at line {next_idx + 1}")
                        item_map.update(extra)
                    items.append(item_map)
                    continue

                items.append(parse_scalar(after_dash))
                idx += 1

            return items, idx

        # mapping
        mapping: Dict[str, Any] = {}
        while True:
            idx = next_nonempty(idx)
            if idx >= n:
                break
            line = lines[idx]
            ind = indent_of(line)
            if ind < current_indent:
                break
            if ind != current_indent:
                raise ParseError(
                    f"{path}: mapping entry has unexpected indent (expected {current_indent}, got {ind})"
                )

            if line.lstrip(" ").startswith("- "):
                break

            if ":" not in line:
                raise ParseError(f"{path}: expected key: value at line {idx + 1}")

            key_raw, rest = line.split(":", 1)
            key = key_raw.strip()
            rest = rest.lstrip()

            if rest in ("|", ">"):
                idx += 1
                scalar, idx = parse_block_scalar(idx, current_indent + 2, rest)
                mapping[key] = scalar
                continue

            if rest == "":
                idx += 1
                child, idx = parse_node(idx, current_indent + 2)
                mapping[key] = child
                continue

            mapping[key] = parse_scalar(rest)
            idx += 1

        return mapping, idx

    try:
        node, _ = parse_node(0, 0)
    except ParseError as e:
        raise SystemExit(
            f"Failed to parse YAML scenario (PyYAML not installed). {e}\n"
            "Install PyYAML for full YAML support, or keep scenarios within the supported subset."
        )

    if not isinstance(node, dict):
        raise SystemExit(f"{path}: expected a YAML mapping at the document root.")
    return node


def replace_vars(command: str, variables: Dict[str, str]) -> str:
    result = command
    for key, value in variables.items():
        token = f"{{{{{key}}}}}"
        if token in result:
            result = result.replace(token, value)
    return result


def run_shell(
    command: str,
    datasite: Optional[str],
    variables: Dict[str, str],
    capture: bool = False,
) -> subprocess.CompletedProcess:
    expanded = replace_vars(command, variables)
    env = os.environ.copy()
    env["WORKSPACE_ROOT"] = str(ROOT)
    env["SCENARIO_DATASITE"] = datasite or ""
    java_home = JAVA_HOME_OVERRIDE
    user_path = env.get("SCENARIO_USER_PATH")

    if datasite:
        client_dir = SANDBOX_ROOT / datasite
        if not client_dir.is_dir():
            raise SystemExit(f"Sandbox datasite missing: {client_dir}")

        data_dir = client_dir
        config_path = client_dir / ".syftbox" / "config.json"
        if not config_path.exists():
            raise SystemExit(f"Missing SyftBox config for {datasite}: {config_path}")

        env.update(
            {
                "HOME": str(client_dir),
                "SYFTBOX_EMAIL": datasite,
                "SYFTBOX_DATA_DIR": str(data_dir),
                "SYFTBOX_CONFIG_PATH": str(config_path),
            }
        )
        if java_home:
            env["JAVA_HOME"] = java_home
            env["JAVA_CMD"] = str(Path(java_home) / "bin" / "java")
        # Relax Nextflow Java version checks to allow newer JVMs (e.g., 25)
        env.setdefault("NXF_DISABLE_JAVA_VERSION_CHECK", "true")
        env.setdefault("NXF_IGNORE_JAVA_VERSION", "true")
        env.setdefault("NXF_OPTS", "-Dnxf.java.check=false")
        if user_path:
            env["PATH"] = user_path

        result = subprocess.run(
            ["bash", "-c", expanded],
            cwd=client_dir,
            env=env,
            text=True,
            capture_output=True,
        )
    else:
        if java_home:
            env["JAVA_HOME"] = java_home
            env["JAVA_CMD"] = str(Path(java_home) / "bin" / "java")
        env.setdefault("NXF_DISABLE_JAVA_VERSION_CHECK", "true")
        env.setdefault("NXF_IGNORE_JAVA_VERSION", "true")
        env.setdefault("NXF_OPTS", "-Dnxf.java.check=false")
        result = subprocess.run(
            expanded,
            cwd=ROOT,
            env=env,
            text=True,
            shell=True,
            capture_output=True,
        )

    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)
    if result.returncode != 0:
        raise SystemExit(f"Command failed (exit {result.returncode}): {expanded}")
    return result


def execute_commands(commands: Any, variables: Dict[str, str]):
    if not commands:
        return
    for cmd in commands:
        run_shell(cmd, None, variables)


def run_step(step: Dict[str, Any], variables: Dict[str, str]):
    name = step.get("name")
    if name:
        print(f"\n=== {name} ===")
    datasite = step.get("datasite")
    command = step.get("run")
    capture = step.get("capture")
    expect = step.get("expect")
    wait_for = step.get("wait_for")
    timeout = step.get("timeout", 30)
    assert_no_encrypted = step.get("assert_no_encrypted")
    assert_encrypted = step.get("assert_encrypted")

    # Handle wait_for syntax
    if wait_for:
        import time
        import glob
        expanded = replace_vars(wait_for, variables)
        wait_for_new = bool(step.get("wait_for_new", False))

        if datasite:
            # Wait for file in datasite directory
            client_dir = SANDBOX_ROOT / datasite
            wait_path = client_dir / expanded
        else:
            wait_path = ROOT / expanded

        initial_matches: set[str] = set()
        if wait_for_new and "*" in str(wait_path):
            initial_matches = set(glob.glob(str(wait_path)))

        for i in range(timeout):
            # Support glob patterns
            if '*' in str(wait_path):
                matches = set(glob.glob(str(wait_path)))
                candidates = sorted(matches - initial_matches) if wait_for_new else sorted(matches)
                if candidates:
                    print(f"✓ Found: {candidates[0]}")
                    return
            elif wait_path.exists():
                print(f"✓ Found: {wait_path}")
                return
            time.sleep(1)

        raise SystemExit(f"Timeout waiting for: {wait_path}")

    # Handle assert_no_encrypted: check that files/directory have no SYC1 headers
    if assert_no_encrypted:
        import glob
        expanded = replace_vars(assert_no_encrypted, variables)

        if datasite:
            client_dir = SANDBOX_ROOT / datasite
            check_path = client_dir / expanded
        else:
            check_path = ROOT / expanded

        # Collect files to check
        if check_path.is_dir():
            files_to_check = list(check_path.rglob("*"))
            files_to_check = [f for f in files_to_check if f.is_file()]
        elif '*' in str(check_path):
            files_to_check = [Path(p) for p in glob.glob(str(check_path))]
        else:
            files_to_check = [check_path] if check_path.exists() else []

        # Check each file for SYC1 header
        encrypted_files = []
        for file_path in files_to_check:
            try:
                with open(file_path, 'rb') as f:
                    header = f.read(4)
                    if header == b'SYC1':
                        encrypted_files.append(file_path)
            except Exception:
                pass  # Skip files we can't read

        if encrypted_files:
            print(f"ERROR: Found {len(encrypted_files)} encrypted file(s):")
            for f in encrypted_files:
                print(f"  - {f}")
            raise SystemExit(f"Encrypted files found in: {check_path}")

        print(f"✓ No encrypted files in: {expanded}")
        return

    # Handle assert_encrypted: check that files have SYC1 headers
    if assert_encrypted:
        import glob
        expanded = replace_vars(assert_encrypted, variables)

        if datasite:
            client_dir = SANDBOX_ROOT / datasite
            check_path = client_dir / expanded
        else:
            check_path = ROOT / expanded

        # Collect files to check
        if check_path.is_dir():
            files_to_check = list(check_path.rglob("*"))
            files_to_check = [f for f in files_to_check if f.is_file()]
        elif '*' in str(check_path):
            files_to_check = [Path(p) for p in glob.glob(str(check_path))]
        else:
            files_to_check = [check_path] if check_path.exists() else []

        if not files_to_check:
            raise SystemExit(f"No files found to check: {check_path}")

        # Check each file for SYC1 header
        unencrypted_files = []
        for file_path in files_to_check:
            try:
                with open(file_path, 'rb') as f:
                    header = f.read(4)
                    if header != b'SYC1':
                        unencrypted_files.append(file_path)
            except Exception:
                unencrypted_files.append(file_path)  # Assume unencrypted if can't read

        if unencrypted_files:
            print(f"ERROR: Found {len(unencrypted_files)} unencrypted file(s):")
            for f in unencrypted_files:
                print(f"  - {f}")
            raise SystemExit(f"Unencrypted files found in: {check_path}")

        print(f"✓ All files encrypted in: {expanded}")
        return

    if not command:
        raise SystemExit("Each step must define a 'run', 'wait_for', 'assert_no_encrypted', or 'assert_encrypted' command.")

    result = run_shell(command, datasite, variables, capture=bool(capture))

    if expect:
        if expect not in result.stdout:
            raise SystemExit(f"Expected substring '{expect}' not found in output.")

    if capture:
        variables[capture] = result.stdout.strip()
        print(f"[captured {capture}]")


def main():
    parser = argparse.ArgumentParser(description="Run YAML scenarios against sandbox datasites.")
    parser.add_argument("scenario", type=Path, help="Path to YAML scenario file")
    args = parser.parse_args()

    scenario = load_scenario(args.scenario)
    variables: Dict[str, str] = {"workspace": str(ROOT)}

    setup = scenario.get("setup", {})
    execute_commands(setup.get("server"), variables)
    execute_commands(setup.get("clients"), variables)

    steps = scenario.get("steps", [])
    if not isinstance(steps, list):
        raise SystemExit("Scenario 'steps' must be a list.")

    for step in steps:
        if not isinstance(step, dict):
            raise SystemExit("Each step must be a mapping.")
        run_step(step, variables)

    print("\nScenario completed successfully.")


if __name__ == "__main__":
    main()
