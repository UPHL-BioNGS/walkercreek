process IRMA_ABRICATE_REPORT {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(irma_tsv), path(abricate_tsv)

    output:
    tuple val(meta), path("*.combined.typing.tsv"), emit: tsv_combined
    path "versions.yml", optional: true, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_tsv = "${prefix}.combined.typing.tsv"

    """
    python - <<'PY'
import csv

irma_file = "${irma_tsv}"
abricate_file = "${abricate_tsv}"
out_file = "${output_tsv}"

# Defaults used when a record is missing from either side
IRMA_DEFAULT_TYPE = "No IRMA type"
IRMA_DEFAULT_SUBTYPE = "No IRMA subtype"
ABR_DEFAULT_TYPE = "No abricate type"
ABR_DEFAULT_SUBTYPE = "No abricate subtype"

def norm(v):
    return (v or "").strip()

# Read IRMA table (primary table: we preserve all of these rows)
irma_rows = {}
with open(irma_file, "r", newline="") as f:
    reader = csv.DictReader(f, delimiter="\\t")
    if not reader.fieldnames or "Sample" not in reader.fieldnames:
        raise ValueError(f"IRMA TSV missing expected header 'Sample': {irma_file}")
    for row in reader:
        sample = norm(row.get("Sample"))
        if not sample:
            continue
        irma_rows[sample] = {
            "IRMA_type": norm(row.get("IRMA_type")) or IRMA_DEFAULT_TYPE,
            "IRMA_subtype": norm(row.get("IRMA_subtype")) or IRMA_DEFAULT_SUBTYPE,
        }

# Read ABRICATE table (optional/secondary)
abricate_rows = {}
with open(abricate_file, "r", newline="") as f:
    reader = csv.DictReader(f, delimiter="\\t")
    if not reader.fieldnames or "Sample" not in reader.fieldnames:
        raise ValueError(f"ABRICATE TSV missing expected header 'Sample': {abricate_file}")
    for row in reader:
        sample = norm(row.get("Sample"))
        if not sample:
            continue
        abricate_rows[sample] = {
            "abricate_InsaFlu_type": norm(row.get("abricate_InsaFlu_type")) or ABR_DEFAULT_TYPE,
            "abricate_InsaFlu_subtype": norm(row.get("abricate_InsaFlu_subtype")) or ABR_DEFAULT_SUBTYPE,
        }

# Write a LEFT JOIN from IRMA -> ABRICATE
fieldnames = [
    "Sample",
    "IRMA_type",
    "IRMA_subtype",
    "abricate_InsaFlu_type",
    "abricate_InsaFlu_subtype",
]

with open(out_file, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\\t")
    writer.writeheader()

    for sample in sorted(irma_rows.keys()):
        ir = irma_rows[sample]
        ab = abricate_rows.get(sample, {
            "abricate_InsaFlu_type": ABR_DEFAULT_TYPE,
            "abricate_InsaFlu_subtype": ABR_DEFAULT_SUBTYPE,
        })

        writer.writerow({
            "Sample": sample,
            "IRMA_type": ir["IRMA_type"],
            "IRMA_subtype": ir["IRMA_subtype"],
            "abricate_InsaFlu_type": ab["abricate_InsaFlu_type"],
            "abricate_InsaFlu_subtype": ab["abricate_InsaFlu_subtype"],
        })
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
