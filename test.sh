#!/usr/bin/env bash

cd $(dirname 0)

CASE=${1:-reads0}

python3 ./testerka.py ./example_files1/reference.fasta "./example_files1/$CASE.fasta" "./example_files1/$CASE.txt" ./app/main.py
