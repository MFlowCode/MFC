#!/usr/bin/env python3
"""
IB accuracy limits: the 2D_ibm_conduction_annulus deck (steady conduction between an isothermal circle and an isothermal STL ring,
exact ln(r) solution), run from this directory so that order_study.py can refine it. Takes the annulus arguments (--cells_per_R,
--interp, --steps); see ../2D_ibm_conduction_annulus/case.py.

    ./mfc.sh run examples/2D_ibm_accuracy_limits/case.py -n 1 -- --cells_per_R 20 --interp 2
"""

import json
import os
import runpy

ANNULUS = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "2D_ibm_conduction_annulus")
case = runpy.run_path(os.path.join(ANNULUS, "case.py"))["case"]
case["stl_models(1)%model_filepath"] = os.path.abspath(os.path.join(ANNULUS, "ring.stl"))

if __name__ == "__main__":
    print(json.dumps(case))
