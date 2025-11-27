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
            raise SystemExit(
                "PyYAML is required for non-JSON scenarios. Install it with 'pip install pyyaml'."
            ) from exc
        return yaml.safe_load(text)  # type: ignore[arg-type]


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

        if datasite:
            # Wait for file in datasite directory
            client_dir = SANDBOX_ROOT / datasite
            wait_path = client_dir / expanded
        else:
            wait_path = ROOT / expanded

        for i in range(timeout):
            # Support glob patterns
            if '*' in str(wait_path):
                matches = glob.glob(str(wait_path))
                if matches:
                    print(f"✓ Found: {matches[0]}")
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
