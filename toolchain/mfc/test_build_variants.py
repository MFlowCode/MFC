"""Tests for the concurrent build variants in .github/workflows/common/build.sh.

The suite needs four kinds of simulation binary (base, two chemistry mechanisms,
state-dependent EOS), and a Frontier job built them one after another: 67-89 minutes with CCE on
GPUs, of a two-hour allocation, so a retry after a late failure could not finish.
They stage under different slugs, the allocation is a whole node, and the AMD GPU
lane has long built them as separate concurrent jobs. Here they build side by side
inside the one allocation (the g1 partition caps a user at 16 nodes, so more jobs
would only queue), and the full build runs afterwards as a sweep: it finds
everything in place, and builds whatever a variant missed.
"""

import os
import shutil
import stat
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# Logs each call with its start and end time; syscheck installs a passing probe.
MFC_STUB = """#!/bin/bash
start=$(date +%s.%N)
for a in "$@"; do
  if [ "$a" = "syscheck" ]; then
    echo "$start $(date +%s.%N) $*" >> trace.log
    mkdir -p build/install/gpu/bin
    printf '#!/bin/bash\\nexit 0\\n' > build/install/gpu/bin/syscheck
    chmod +x build/install/gpu/bin/syscheck
    exit 0
  fi
done
echo "output of: $*"
case " $* " in *" $FAIL_ARG "*) rc=1 ;; *) rc=0 ;; esac
case " $* " in *" -o "*|*" build -j "*) sleep 1 ;; esac
echo "$start $(date +%s.%N) $*" >> trace.log
exit $rc
"""


def _run(tmp_path, cluster="frontier", variant="", fail_arg="none"):
    shutil.copytree(REPO / ".github", tmp_path / ".github")
    binz = tmp_path / "bin"
    binz.mkdir()
    passthrough = '#!/bin/bash\nwhile [ "${1:0:1}" = "-" ]; do shift; case "$1" in [0-9]*) shift;; esac; done\nexec "$@"\n'
    for name, text in (("mpirun", passthrough), ("srun", passthrough)):
        (binz / name).write_text(text)
        (binz / name).chmod(0o755)
    (tmp_path / "mfc.sh").write_text(MFC_STUB)
    (tmp_path / "mfc.sh").chmod((tmp_path / "mfc.sh").stat().st_mode | stat.S_IEXEC)
    env = {
        **os.environ,
        "PATH": f"{binz}:{os.environ['PATH']}",
        "MFC_CI_STATE_DIR": str(tmp_path / "state"),
        "MFC_BUILD_RETRY_DELAY": "0",
        "job_device": "gpu",
        "job_interface": "acc",
        "job_shard": "",
        "job_cluster": cluster,
        "job_variant": variant,
        "SLURMD_NODENAME": "frontier1234",
        "SLURM_JOB_ID": "123456",
        "FAIL_ARG": fail_arg,
    }
    result = subprocess.run(["bash", ".github/workflows/common/build.sh"], capture_output=True, text=True, cwd=tmp_path, env=env, check=False, timeout=120)
    calls = [(float(s), float(e), cmd) for s, e, cmd in (line.split(" ", 2) for line in (tmp_path / "trace.log").read_text().splitlines())]
    return result, sorted(calls)


def _is(kind, cmd):
    return {
        "base": cmd.startswith("build -j"),
        "chem": "-o Chemistry" in cmd and ".yaml" not in cmd,
        "chem-h2o2": "-o Chemistry h2o2.yaml" in cmd,
        "chem-carbon": "-o Chemistry carbon_gasphase_reduced_gri11.yaml" in cmd,
        "eos": "-o eos=mie_gruneisen" in cmd,
        "sweep": cmd.startswith("test") and " -o " not in cmd,
    }[kind]


SIDE_BY_SIDE = ("base", "chem-h2o2", "chem-carbon", "eos")


def test_the_variants_build_side_by_side_then_the_sweep(tmp_path):
    result, calls = _run(tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    variants = [c for c in calls if any(_is(k, c[2]) for k in SIDE_BY_SIDE)]
    assert sorted(k for c in variants for k in SIDE_BY_SIDE if _is(k, c[2])) == sorted(SIDE_BY_SIDE), calls
    assert max(s for s, _, _ in variants) < min(e for _, e, _ in variants), "the variants did not overlap"
    sweep = [c for c in calls if _is("sweep", c[2])]
    assert len(sweep) == 1 and sweep[0][0] >= max(e for _, e, _ in variants), "the full build must run once, after the variants"
    assert "syscheck" in calls[0][2], "the probe is still built first"


def test_each_variants_output_is_labelled_and_kept_in_its_own_log(tmp_path):
    result, _ = _run(tmp_path)
    for name, arg in (("base", "build -j 8 --gpu acc"), ("chem-h2o2", "-o Chemistry h2o2.yaml"), ("chem-carbon", "-o Chemistry carbon_"), ("eos", "-o eos=mie_gruneisen")):
        assert any(line.startswith(f"[{name}] ") and arg in line for line in result.stdout.splitlines()), f"no labelled output for {name}"
        assert arg in (tmp_path / f"build-gpu-acc-{name}.log").read_text()


def test_a_failed_variant_fails_the_job_after_the_others_finish(tmp_path):
    result, calls = _run(tmp_path, fail_arg="h2o2.yaml")
    assert result.returncode == 1
    assert not [c for c in calls if _is("sweep", c[2])], "the sweep must not run after a failed variant"
    assert all([c for c in calls if _is(k, c[2])] for k in ("base", "chem-carbon", "eos")), "the other variants must be allowed to finish"
    assert len([c for c in calls if _is("chem-h2o2", c[2])]) == 2, "the failed variant gets its one retry"


def test_a_job_given_one_variant_builds_only_that_variant(tmp_path):
    """The AMD GPU lane submits one job per variant (test.yml); those jobs are unchanged."""
    result, calls = _run(tmp_path, cluster="frontier_amd", variant="chem")
    assert result.returncode == 0, result.stdout + result.stderr
    builds = [c[2] for c in calls if "syscheck" not in c[2]]
    assert len(builds) == 1 and _is("chem", builds[0]), builds
