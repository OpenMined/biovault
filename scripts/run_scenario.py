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
SBENV = ROOT / "sbenv" / "sbenv"


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
    java_home = env.get("SCENARIO_JAVA_HOME")
    user_path = env.get("SCENARIO_USER_PATH")

    if datasite:
        client_dir = SANDBOX_ROOT / datasite
        if not client_dir.is_dir():
            raise SystemExit(f"Sandbox datasite missing: {client_dir}")
        shell_parts = [
            'set -eo pipefail',
            f'eval "$({SBENV} activate --quiet)"',
        ]
        if java_home:
            shell_parts.append(f'export JAVA_HOME="{java_home}"')
        if user_path:
            shell_parts.append(f'export PATH="{user_path}"')
        shell_parts.append(expanded)
        shell = '; '.join(shell_parts)
        result = subprocess.run(
            ["bash", "-lc", shell],
            cwd=client_dir,
            env=env,
            text=True,
            capture_output=True,
        )
    else:
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

    if not command:
        raise SystemExit("Each step must define a 'run' or 'wait_for' command.")

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
