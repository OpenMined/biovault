#!/bin/bash

SUBMISSION_YAML=$1
RUN_TEST=${2:-false}
echo "Submission YAML: $SUBMISSION_YAML"
if [ ! -f "$SUBMISSION_YAML" ]; then
    echo "Error: File '$SUBMISSION_YAML' does not exist."
    exit 1
fi

WORKFLOW_FILE=$(yq e '.workflow' "$SUBMISSION_YAML")
WORKFLOW_FILE="$(dirname "$SUBMISSION_YAML")/$WORKFLOW_FILE"

if [ "$RUN_TEST" = true ]; then
    PATIENTS=("TEST")
else
    PATIENTS=($(yq e '.patients[]' "$SUBMISSION_YAML"))
fi

echo "List of patients:"
echo "${PATIENTS[@]}"

ASSETS_DIR="./submission/assets"
mkdir -p "$ASSETS_DIR"

for PATIENT in "${PATIENTS[@]}"; do
    get_patient_data() {
        local patient=$1
        local patient_data=$(yq e ".patient.$patient" pt.yaml)
        local ref_version=$(echo "$patient_data" | yq e '.ref_version' -)
        local ref=$(echo "$patient_data" | yq e '.ref' -)
        local ref_index=$(echo "$patient_data" | yq e '.ref_index' -)
        local aligned=$(echo "$patient_data" | yq e '.aligned' -)
        local aligned_index=$(echo "$patient_data" | yq e '.aligned_index' -)
        
        echo "$ref_version" "$ref" "$ref_index" "$aligned" "$aligned_index"
    }

    read REF_VERSION REF REF_INDEX ALIGNED ALIGNED_INDEX <<< $(get_patient_data "$PATIENT")

    RESULTS_DIR="./submission/results/$PATIENT"
    mkdir -p "$RESULTS_DIR"
    echo "Processing patient: $PATIENT"
    echo "nextflow run template.nf --ref_version $REF_VERSION --ref $REF --ref_index $REF_INDEX --aligned $ALIGNED --aligned_index $ALIGNED_INDEX --assets_dir $ASSETS_DIR --results_dir $RESULTS_DIR --work_flow_file $WORKFLOW_FILE -with-docker"
    nextflow run template.nf \
        --patient_id $PATIENT \
        --ref_version $REF_VERSION \
        --ref $REF \
        --ref_index $REF_INDEX \
        --aligned $ALIGNED \
        --aligned_index $ALIGNED_INDEX \
        --assets_dir $ASSETS_DIR \
        --results_dir $RESULTS_DIR \
        --work_flow_file $WORKFLOW_FILE \
        -with-docker

done
