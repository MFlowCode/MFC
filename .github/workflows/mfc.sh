#!/bin/bash
echo "Mock benchmark running..."
mkdir -p bench
echo "name: test_benchmark" > bench/test.yaml
echo "time: 1.23" >> bench/test.yaml
echo "Mock benchmark completed"
