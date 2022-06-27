import copy

import common

from mfc.case import *

from tests.case import TestCase, TCV, TCVS
from tests.case import create_case
from tests.case import CaseGeneratorStack


def get_dimensions():
    r = []

    for dimInfo in [
        (["x"], TCVS(
            trace="1D (m=299,n=0,p=0)",
            variations=[
                TCV("domain.cells.x", 299),
            ]
        ), PatchGeometry.D1_LINE_SEGMENT),
        (["x", "y"], TCVS(
            trace="2D (m=49,n=39,p=0)",
            variations=[
                TCV("domain.cells.x", 49),
                TCV("domain.cells.y", 39),
            ]
        ), PatchGeometry.D2_RECTANGLE),
        (["x", "y", "z"], TCVS(
            trace="3D (m=24,n=24,p=24)",
            variations=[
                TCV("domain.cells.x", 24),
                TCV("domain.cells.y", 24),
                TCV("domain.cells.z", 24),
            ]
        ), PatchGeometry.D3_CUBOID)
    ]:
        dimParams: TCVS = copy.deepcopy(dimInfo[1])

        for dimCmp in dimInfo[0]:
            dimParams.adds([
                TCV(f"domain.domain.{dimCmp}.begin", 0.E+00),
                TCV(f"domain.domain.{dimCmp}.end",   1.E+00)
            ])
        
        if "z" in dimInfo[0]:
            dimParams.adds([
                TCV(f"patches[0].centroid.z", 0.05),
                TCV(f"patches[0].length.z",   0.1),
                TCV(f"patches[1].centroid.z", 0.45),
                TCV(f"patches[1].length.z",   0.7),
                TCV(f"patches[2].centroid.z", 0.9),
                TCV(f"patches[2].length.z",   0.2),
            ])
        elif "y" in dimInfo[0]:
            dimParams.adds([
                TCV(f"patches[0].centroid.y", 0.05),
                TCV(f"patches[0].length.y",   0.1),
                TCV(f"patches[1].centroid.y", 0.45),
                TCV(f"patches[1].length.y",   0.7),
                TCV(f"patches[2].centroid.y", 0.9),
                TCV(f"patches[2].length.y",   0.2)
            ])
        else:
            dimParams.adds([
                TCV(f"patches[0].centroid.x", 0.05),
                TCV(f"patches[0].length.x",   0.1),
                TCV(f"patches[1].centroid.x", 0.45),
                TCV(f"patches[1].length.x",   0.7),
                TCV(f"patches[2].centroid.x", 0.9),
                TCV(f"patches[2].length.x",   0.2)
            ])

        for patchID in range(3):
            dimParams.add(TCV(f"patches[{patchID}].geometry", dimInfo[2]))

            if "z" in dimInfo[0]:
                dimParams.adds([
                    TCV(f"patches[{patchID}].centroid.y", 0.5),
                    TCV(f"patches[{patchID}].length.y",   1),
                    TCV(f"patches[{patchID}].centroid.x", 0.5),
                    TCV(f"patches[{patchID}].length.x",   1)
                ])

            elif "y" in dimInfo[0]:
                dimParams.adds([
                    TCV(f"patches[{patchID}].centroid.x", 0.5),
                    TCV(f"patches[{patchID}].length.x",   1)
                ])

            for dimCmp in dimInfo[0]:
                dimParams.add(TCV(f"patches[{patchID}].velocity.{dimCmp}", 0.0))
            
        r.append((dimInfo[0], dimParams))

    return r


def generate_cases() -> list:
    stack, cases = CaseGeneratorStack(), []

    for dims, dimParams in get_dimensions():
        stack.push(dimParams)

        for bc in [ -1, -2, -4, -5, -6, -7, -8, -9, -10, -11, -12, -3 ]:
            bc_variations = []
            for dim in dims:
                bc_variations = bc_variations + [
                    TCV(f"algorithm.boundary.{dim}.begin", bc),
                    TCV(f"algorithm.boundary.{dim}.end",   bc),                    
                 ]

            stack.push(TCVS(f"bc={bc}", bc_variations))

            cases.append(create_case(stack))

            if bc != -3: # Use bc = 3 henceforth
                stack.pop()

        for weno_order in [3, 5]:
            stack.push(TCVS(f"weno_order={weno_order}", [
                TCV("algorithm.weno.order", weno_order)
            ]))

            for mapped_weno, mp_weno in [('F', 'F'), ('T', 'F'), ('F', 'T')]:
                stack.push(TCVS(f"(mapped_weno={mapped_weno},mp_weno={mp_weno})", [
                    TCV("algorithm.weno.mapped", mapped_weno == 'T'),
                    TCV("algorithm.weno.mp",     mp_weno == 'T')
                ]))

                if not (mp_weno == 'T' and weno_order != 5):
                    cases.append(create_case(stack))

                stack.pop()

            stack.pop()

        for num_fluids in [1, 2]:
            stack.push(TCVS(f"num_fluids={num_fluids}", [TCV("fluids.count", num_fluids)]))

            if num_fluids == 2:
                stack.push(TCVS(None, [
                    # 2nd Fluid
                    TCV("fluids[1].gamma",  2.5),
                    TCV("fluids[1].pi_inf", 0.0),
                    # Patches
                    TCV("patches[0].alpha_rho",    [0.81, 0.19  ]),
                    TCV("patches[0].alpha",        [0.9,  0.1   ]),
                    TCV("patches[1].alpha_rho",    [0.25, 0.25  ]),
                    TCV("patches[1].alpha",        [0.5,  0.5   ]),
                    TCV("patches[2].alpha_rho",    [0.08, 0.0225]),
                    TCV("patches[2].alpha",        [0.2,  0.8   ]),
                ]))

            for riemann_solver in [RiemannSolver.HLL, RiemannSolver.HLLC]:
                stack.push(TCVS(f"riemann_solver={riemann_solver.value}", [
                    TCV("algorithm.riemann_solver", riemann_solver)
                ]))

                cases.append(create_case(stack, TCVS("mixture_err=T", [
                    TCV("algorithm.mixture_err", True)
                ])))
                
                cases.append(create_case(stack, TCVS("avg_state=1", [
                    TCV("algorithm.avg_state", AverageStateEvaluation.ROE_MEAN)
                ])))

                cases.append(create_case(stack, TCVS("wave_speeds=2", [
                    TCV("algorithm.wave_speeds", WaveSpeedEstimation.PRESSURE_VELOCITY)
                ])))

                if num_fluids == 2:
                    if riemann_solver == RiemannSolver.HLLC:
                        cases.append(create_case(stack, TCVS("alt_soundspeed=T", [
                            TCV("algorithm.alt_soundspeed", True)
                        ])))

                    cases.append(create_case(stack, TCVS("mpp_lim=T", [
                        TCV("algorithm.mpp_lim", True)
                    ])))

                stack.pop()

            if num_fluids == 2:
                stack.pop()

            stack.pop()

        if len(dims) == 3:
            cases.append(create_case(stack, TCVS(f"ppn=2,m=29,n=29,p=49", [
                TCV("domain.cells", f"Cells(29, 29, 49)")
            ]), ppn=2))
        else:
            cases.append(create_case(stack, TCVS(f"ppn=2", []), ppn=2))

        stack.push(TCVS("", [
            TCV("domain.time.dt", str([1e-07, 1e-06, 1e-06][len(dims)-1]))
        ]))

        stack.push(TCVS(f"bubbles=T", [
            TCV("bubbles.bubbles", True)
        ]))

        stack.push(TCVS("", [
            TCV("bubbles.number", 3),
            # Fluids
            TCV("fluids[0].gamma",   0.16     ), 
            TCV("fluids[0].pi_inf",  3515.0   ),
            TCV("fluids[0].mul0",    0.001002 ),
            TCV("fluids[0].ss",      0.07275  ),
            TCV("fluids[0].pv",      2338.8   ),
            TCV("fluids[0].gamma_v", 1.33     ),
            TCV("fluids[0].M_v",     18.02    ),
            TCV("fluids[0].mu_v",    8.816e-06),
            TCV("fluids[0].k_v",     0.019426 ),
            TCV("fluids[1].gamma",         2.5), 
            TCV("fluids[1].pi_inf",        0.0), 
            TCV("fluids[1].gamma_v",       1.4),
            TCV("fluids[1].M_v",         28.97),
            TCV("fluids[1].mu_v",      1.8e-05),
            TCV("fluids[1].k_v",       0.02556), 
            # Patches
            TCV("patches[0].alpha_rho", [0.999999999999]), 
            TCV("patches[0].alpha",     [1e-12]),
            TCV("patches[0].pressure",  1.0),
            TCV("patches[1].alpha_rho", [0.96]), 
            TCV("patches[1].alpha",     [4e-02]),  
            TCV("patches[1].pressure",  1.0),
            TCV("patches[2].alpha_rho", [0.999999999999]),
            TCV("patches[2].alpha",     [1e-12]), 
            TCV("patches[2].pressure",  1.0)
        ]))

        stack.push(TCVS("Monopole=T", [
            TCV("acoustic.monopole", True)
        ]))

        # TODO: mono(1)%loc(1) was never set
        if len(dims) >= 2:
            stack.push(TCVS("", [
                TCV("acoustic.monopoles[0].location.y", 0.5)
            ]))

        if len(dims) >= 3:
            stack.push(TCVS("", [
                TCV("acoustic.monopoles[0].location.z", 0.5),
                TCV("acoustic.monopoles[0].support", "AcousticSpacialSupport.D3_FINITE_LINE_PATCH")
            ]))

        for polytropic in ['T', 'F']:
            stack.push(TCVS(f"polytropic={polytropic}", [
                TCV("bubbles.polytropic", polytropic == 'T')
            ]))

            for bubble_model in [
                BubbleModel.RAYLEIGH_PLESSET,
                BubbleModel.KELLER_MIKSIS
            ]:
                stack.push(TCVS(f"bubble_model={bubble_model.value}", [
                    TCV("bubbles.model", bubble_model)
                ]))

                if not (polytropic == 'F' and bubble_model == BubbleModel.RAYLEIGH_PLESSET):
                    cases.append(create_case(stack))

                stack.pop()

            stack.pop()

        stack.push(TCVS("", [
            TCV("bubbles.polytropic", True),
            TCV("bubbles.model",      BubbleModel.KELLER_MIKSIS)
        ]))

        cases.append(create_case(stack, TCVS("nb=1", [
            TCV("bubbles.number", 1)
        ])))

        stack.push(TCVS("qbmm=T", [
            TCV("bubbles.qbmm", True)
        ]))
        cases.append(create_case(stack))
        
        stack.push(TCVS("bubble_model=3", [
            TCV("bubbles.model", BubbleModel.RAYLEIGH_PLESSET)
        ]))
        cases.append(create_case(stack))

        for i in range(6):
            stack.pop()

        if len(dims) >= 2:
            stack.pop()

        if len(dims) >= 3:
            stack.pop()

        for i in range(3):
            stack.pop()

    # Sanity Check 1

    if stack.size() != 0:
        raise common.MFCException("generate_cases: stack isn't fully pop'ed")

    # Sanity Check 2
    uuids = [ case.get_uuid() for case in cases ]
    l1, l2 = len(uuids), len(set(uuids))
    if l1 != l2:
        raise common.MFCException(f"generate_cases: uuids aren't unique ({l1} cases but {l2} unique uuids)")

    return cases
