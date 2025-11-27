#!/usr/bin/python3
import os
import sys

def main(input_files, output_file):
    with open(output_file, "w") as out:
        for path in input_files:
            abs_path = os.path.abspath(path)
            out.write(abs_path + "\n")

if __name__ == "__main__":
    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    main(input_files, output_file)