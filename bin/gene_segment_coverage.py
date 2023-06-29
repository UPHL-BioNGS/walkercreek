import os

root_directory = 'workDir/irma/'

output_file = 'gene_segment_coverage.tsv'

# Recursive function to find all *-coverage.txt files
def find_coverage_files(directory):
    coverage_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('-coverage.txt'):
                coverage_files.append(os.path.join(root, file))
    return coverage_files

# Calculate average coverage for each Reference_Name
def calculate_average_coverage(coverage_file):
    reference_name = os.path.basename(coverage_file).split('-')[0]
    total_rows = 0
    total_coverage = 0.0

    with open(coverage_file, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            data = line.strip().split('\t')
            coverage = float(data[7])
            total_coverage += coverage
            total_rows += 1

    average_coverage = total_coverage / total_rows
    return reference_name, average_coverage

# Find coverage files and calculate average coverage
coverage_files = find_coverage_files(os.path.join(root_directory, 'meta.id'))
results = {}

for coverage_file in coverage_files:
    reference_name, average_coverage = calculate_average_coverage(coverage_file)
    if reference_name not in results:
        results[reference_name] = []
    results[reference_name].append(average_coverage)

# Write results to the output file
with open(output_file, 'w') as file:
    file.write('Reference_Name\tAverage_Coverage\n')
    for reference_name, coverages in results.items():
        average_coverage = sum(coverages) / len(coverages)
        file.write(f'{reference_name}\t{average_coverage}\n')

print(f"Average coverage data written to {output_file}")



