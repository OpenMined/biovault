// BioVault multiparty workflow v0.1.0
// Simple 3-party collaborative sum for testing multiparty UX

nextflow.enable.dsl=2

workflow USER {
    take:
        context
        role        // contributor1, contributor2, or aggregator
        session_id  // unique session identifier
        step_id     // which step to execute
        input_dir   // directory with inputs from other parties (optional)

    main:
        if (step_id == "generate") {
            def result = generate_numbers(role, session_id)
            emit:
                numbers = result.numbers
        } else if (step_id == "aggregate") {
            def result = aggregate_sum(session_id, input_dir)
            emit:
                result = result.sum_result
        }
}

// Generate random numbers (for contributors)
process generate_numbers {
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        val role
        val session_id

    output:
        path "numbers.json", emit: numbers

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Generate 5 random numbers between 1 and 100
    python3 << 'PYEOF'
import json
import random

numbers = [random.randint(1, 100) for _ in range(5)]
result = {
    "role": "${role}",
    "session_id": "${session_id}",
    "numbers": numbers,
    "sum": sum(numbers)
}

with open("numbers.json", "w") as f:
    json.dump(result, f, indent=2)

print(f"Generated numbers: {numbers} (sum: {sum(numbers)})")
PYEOF
    """
}

// Aggregate contributions from all parties
process aggregate_sum {
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        val session_id
        path input_dir

    output:
        path "result.json", emit: sum_result

    script:
    """
    #!/bin/bash
    set -euo pipefail

    python3 << 'PYEOF'
import json
import os
from pathlib import Path

input_path = Path("${input_dir}")
all_numbers = []
contributions = []

# Read all contribution files
for json_file in input_path.glob("*.json"):
    try:
        with open(json_file) as f:
            data = json.load(f)
            if "numbers" in data:
                all_numbers.extend(data["numbers"])
                contributions.append({
                    "role": data.get("role", "unknown"),
                    "numbers": data["numbers"],
                    "partial_sum": data.get("sum", sum(data["numbers"]))
                })
    except Exception as e:
        print(f"Warning: Could not read {json_file}: {e}")

total_sum = sum(all_numbers)

result = {
    "session_id": "${session_id}",
    "contributions": contributions,
    "all_numbers": all_numbers,
    "total_sum": total_sum,
    "count": len(all_numbers)
}

with open("result.json", "w") as f:
    json.dump(result, f, indent=2)

print(f"Aggregated {len(contributions)} contributions")
print(f"Total numbers: {len(all_numbers)}")
print(f"Total sum: {total_sum}")
PYEOF
    """
}
