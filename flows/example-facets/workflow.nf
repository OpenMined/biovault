// Example facet-aware BioVault workflow.
//
// BioVault appends participant facets as extra samplesheet columns. For example:
// participant_id,genotype_file,country,split
// 100001,/path/100001.txt,Bahamas,island
//
// This workflow keeps those extra fields on each record and aggregates a simple
// variant-count proxy by every facet column. The count here is intentionally
// simple for testing: non-empty, non-comment genotype file lines.

nextflow.enable.dsl=2

workflow USER {
    take:
        context
        participants

    main:
        def participant_work_items = participants.map { record ->
            def facets = record.facets instanceof Map ? record.facets : [:]
            if (facets.isEmpty()) {
                def reserved = ['participant_id', 'genotype_file', 'genotype_path', 'validation', 'grch_build', 'source', 'facets'] as Set
                record.each { key, value ->
                    if (!(key in reserved) && value != null && value.toString().trim()) {
                        facets[key] = value.toString()
                    }
                }
            }
            tuple(
                record.participant_id.toString(),
                file(record.genotype_file),
                facets
            )
        }

        def participant_rows = count_variants(participant_work_items)
        def aggregate_rows = aggregate_facets(participant_rows.collect())

    emit:
        participant_counts = aggregate_rows.participant_counts
        facet_counts = aggregate_rows.facet_counts
}

process count_variants {
    tag { participant_id }

    input:
        tuple val(participant_id), path(genotype_file), val(facets)

    output:
        path "${participant_id}.participant.tsv"

    script:
    def facets_json = groovy.json.JsonOutput.toJson(facets)
    """
    python3 - <<'PY'
    import json
    import pathlib

    participant_id = ${groovy.json.JsonOutput.toJson(participant_id)}
    genotype_path = pathlib.Path(${groovy.json.JsonOutput.toJson(genotype_file.getName())})
    facets = json.loads(${groovy.json.JsonOutput.toJson(facets_json)})

    variant_count = 0
    with genotype_path.open('r', encoding='utf-8', errors='ignore') as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith('#'):
                variant_count += 1

    facet_json = json.dumps(facets, sort_keys=True)
    with open(f"{participant_id}.participant.tsv", "w", encoding="utf-8") as out:
        out.write("participant_id\\tvariant_count\\tfacets_json\\n")
        out.write(f"{participant_id}\\t{variant_count}\\t{facet_json}\\n")
    PY
    """
}

process aggregate_facets {
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path participant_files

    output:
        path "participant_counts.tsv", emit: participant_counts
        path "facet_counts.tsv", emit: facet_counts

    script:
    """
    python3 - <<'PY'
    import collections
    import csv
    import glob
    import json

    participant_rows = []
    facet_names = set()

    for path in sorted(glob.glob("*.participant.tsv")):
        with open(path, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\\t")
            for row in reader:
                facets = json.loads(row.get("facets_json") or "{}")
                facet_names.update(facets.keys())
                participant_rows.append({
                    "participant_id": row["participant_id"],
                    "variant_count": int(row["variant_count"]),
                    "facets": facets,
                })

    facet_names = sorted(facet_names)
    with open("participant_counts.tsv", "w", newline="", encoding="utf-8") as out:
        fieldnames = ["participant_id", "variant_count"] + facet_names
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\\t")
        writer.writeheader()
        for row in participant_rows:
            writer.writerow({
                "participant_id": row["participant_id"],
                "variant_count": row["variant_count"],
                **{name: row["facets"].get(name, "") for name in facet_names},
            })

    grouped = collections.defaultdict(lambda: {"participants": 0, "variant_count": 0})
    for row in participant_rows:
        for name, value in row["facets"].items():
            key = (name, value)
            grouped[key]["participants"] += 1
            grouped[key]["variant_count"] += row["variant_count"]

    with open("facet_counts.tsv", "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\\t")
        writer.writerow(["facet_name", "facet_value", "participants", "variant_count"])
        for (name, value), counts in sorted(grouped.items()):
            writer.writerow([name, value, counts["participants"], counts["variant_count"]])
    PY
    """
}
