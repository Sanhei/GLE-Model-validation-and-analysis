#!/bin/bash
head -10000000 0.001000.txt > test.txt
awk 'NR % 10 == 0' test.txt > test2.txt
