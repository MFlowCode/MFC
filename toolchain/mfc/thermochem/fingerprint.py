"""Content identities for mechanisms and the MFC-owned generator."""

import hashlib
import json
from functools import lru_cache
from pathlib import Path


def mechanism_fingerprint(solution):
    """Hash the phase, species and reactions independently of the current state."""
    import cantera as ct

    data = dict(solution.input_data)
    data.pop("state", None)
    data["species"] = [species.input_data for species in solution.species()]
    data["reactions"] = [reaction.input_data for reaction in solution.reactions()]
    data["cantera-version"] = ct.__version__
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()


@lru_cache(maxsize=1)
def generator_fingerprint():
    """Keep builds from different generator revisions in separate staging trees."""
    digest = hashlib.sha256()
    root = Path(__file__).parent
    for name in ("__init__.py", "fortran.py", "expressions.py", "module.fpp.mako"):
        digest.update((root / name).read_bytes())
    return digest.hexdigest()
