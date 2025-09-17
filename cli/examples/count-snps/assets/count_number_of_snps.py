#!/usr/bin/env python3
import sys


def count_snps(input_stream):
    """
    Count the number of SNPs in the input stream.
    Assumes each line represents a SNP entry.
    """
    count = 0
    for line in input_stream:
        line = line.strip()
        if line and not line.startswith("#"):  # Skip empty lines and comments
            count += 1
    return count


if __name__ == "__main__":
    # Read from stdin
    count = count_snps(sys.stdin)
    print(count, end="")
