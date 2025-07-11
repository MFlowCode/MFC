# Continuouscl Benchmarking Results

Generated on: 2025-07-10 00:46:18 UTC

## System Information
- Host: github-runner
- CPUs: 2
- MHz per CPU: 2000

## Benchmark Results

| Test Case | Metric | Time (seconds) | Time (nanoseconds) |
|-----------|--------|----------------|-------------------|
| 5eq_rk3_weno3_hllc | simulation_time | 553.746147 | 553746146600 |
| 5eq_rk3_weno3_hllc | grind_time | 0.000000 | 311 |
| hypo_hll | simulation_time | 1158.006229 | 1158006229400 |
| hypo_hll | grind_time | 0.000000 | 373 |
| ibm | simulation_time | 1560.666652 | 1560666652400 |
| ibm | grind_time | 0.000001 | 1329 |
| viscous_weno5_sgb_acoustic | simulation_time | 968.760392 | 968760392000 |
| viscous_weno5_sgb_acoustic | grind_time | 0.000001 | 1052 |

## Raw Data

```json
{
  "context": {
    "date": "2025-07-10T00:46:18.506114",
    "host_name": "github-runner",
    "executable": "mfc_benchmark",
    "num_cpus": 2,
    "mhz_per_cpu": 2000,
    "cpu_scaling_enabled": false,
    "caches": []
  },
  "benchmarks": [
    {
      "name": "5eq_rk3_weno3_hllc/simulation_time",
      "family_index": 0,
      "per_family_instance_index": 0,
      "run_name": "5eq_rk3_weno3_hllc/simulation_time",
      "run_type": "iteration",
      "repetitions": 1,
      "repetition_index": 0,
      "threads": 1,
      "iterations": 1,
      "real_time": 553746146600.0,
      "cpu_time": 553746146600.0,
      "time_unit": "ns"
    },
    {
      "name": "5eq_rk3_weno3_hllc/grind_time",
      "family_index": 1,
      "per_family_instance_index": 0,
      "run_name": "5eq_rk3_weno3_hllc/grind_time",
      "run_type": "iteration",
      "repetitions": 1,
      "repetition_index": 0,
      "threads": 1,
      "iterations": 1,
      "real_time": 311.46590698,
      "cpu_time": 311.46590698,
      "time_unit": "ns"
    },
    {
      "name": "hypo_hll/simulation_time",
      "family_index": 2,
      "per_family_instance_index": 0,
      "run_name": "hypo_hll/simulation_time",
      "run_type": "iteration",
      "repetitions": 1,
      "repetition_index": 0,
      "threads": 1,
      "iterations": 1,
      "real_time": 1158006229400.0,
      "cpu_time": 1158006229400.0,
      "time_unit": "ns"
    },
    {
      "name": "hypo_hll/grind_time",
      "family_index": 3,
      "per_family_instance_index": 0,
      "run_name": "hypo_hll/grind_time",
      "run_type": "iteration",
      "repetitions": 1,
      "repetition_index": 0,
      "threads": 1,
      "iterations": 1,
      "real_time": 373.3648501,
      "cpu_time": 373.3648501,
      "time_unit": "ns"
    },
    {
      "name": "ibm/simulation_time",
      "family_index": 4,
      "per_family_instance_index": 0,
      "run_name": "ibm/simulation_time",
      "run_type": "iteration",
      "repetitions": 1,
      "repetition_index": 0,
      "threads": 1,
      "iterations": 1,
      "real_time": 1560666652400.0,
      "cpu_time": 1560666652400.0,
      "time_unit": "ns"
    },
    {
      "name": "ibm/grind_time",
      "family_index": 5,
      "per_family_instance_index": 0,
      "run_name": "ibm/grind_time",
      "run_type": "iteration",
      "repetitions": 1,
      "repetition_index": 0,
      "threads": 1,
      "iterations": 1,
      "real_time": 1329.41865639,
      "cpu_time": 1329.41865639,
      "time_unit": "ns"
    },
    {
      "name": "viscous_weno5_sgb_acoustic/simulation_time",
      "family_index": 6,
      "per_family_instance_index": 0,
      "run_name": "viscous_weno5_sgb_acoustic/simulation_time",
      "run_type": "iteration",
      "repetitions": 1,
      "repetition_index": 0,
      "threads": 1,
      "iterations": 1,
      "real_time": 968760392000.0,
      "cpu_time": 968760392000.0,
      "time_unit": "ns"
    },
    {
      "name": "viscous_weno5_sgb_acoustic/grind_time",
      "family_index": 7,
      "per_family_instance_index": 0,
      "run_name": "viscous_weno5_sgb_acoustic/grind_time",
      "run_type": "iteration",
      "repetitions": 1,
      "repetition_index": 0,
      "threads": 1,
      "iterations": 1,
      "real_time": 1051.98369313,
      "cpu_time": 1051.98369313,
      "time_unit": "ns"
    }
  ]
}
```

---
*Last updated: 2025-07-10T00:46:18.566505*