import typing

from mfc   import common
from .case import TestCase, create_case, CaseGeneratorStack

def get_bc_mods(bc: int, dimInfo):
    params = {}
    for dimCmp in dimInfo[0]:
        params.update({f'bc_{dimCmp}%beg': bc, f'bc_{dimCmp}%end': bc})

    return params


def get_dimensions():
    r = []

    for dimInfo in [ (["x"],           {'m': 299, 'n': 0,  'p': 0},  {"geometry": 1}),
                     (["x", "y"],      {'m': 49,  'n': 39, 'p': 0},  {"geometry": 3}),
                     (["x", "y", "z"], {'m': 24,  'n': 24, 'p': 24}, {"geometry": 9}) ]:
        dimParams = {**dimInfo[1]}

        for dimCmp in dimInfo[0]:
            dimParams.update({
                f"{dimCmp}_domain%beg": 0.E+00, f"{dimCmp}_domain%end": 1.E+00
            })

        dimParams.update(get_bc_mods(-3, dimInfo))

        for patchID in range(1, 3+1):
            dimParams[f"patch_icpp({patchID})%geometry"] = dimInfo[2].get("geometry")

            if "z" in dimInfo[0]:
                dimParams.update({
                    f"patch_icpp({1})%z_centroid":       0.05, f"patch_icpp({1})%length_z":         0.1,
                    f"patch_icpp({2})%z_centroid":       0.45, f"patch_icpp({2})%length_z":         0.7,
                    f"patch_icpp({3})%z_centroid":       0.9,  f"patch_icpp({3})%length_z":         0.2,
                    f"patch_icpp({patchID})%y_centroid": 0.5,  f"patch_icpp({patchID})%length_y":   1,
                    f"patch_icpp({patchID})%x_centroid": 0.5,  f"patch_icpp({patchID})%length_x":   1
                })

            elif "y" in dimInfo[0]:
                dimParams.update({
                    f"patch_icpp({1})%y_centroid":       0.05, f"patch_icpp({1})%length_y":       0.1,
                    f"patch_icpp({2})%y_centroid":       0.45, f"patch_icpp({2})%length_y":       0.7,
                    f"patch_icpp({3})%y_centroid":       0.9,  f"patch_icpp({3})%length_y":       0.2,
                    f"patch_icpp({patchID})%x_centroid": 0.5,  f"patch_icpp({patchID})%length_x": 1
                })
            else:
                dimParams.update({
                    f"patch_icpp({1})%x_centroid": 0.05, f"patch_icpp({1})%length_x": 0.1,
                    f"patch_icpp({2})%x_centroid": 0.45, f"patch_icpp({2})%length_x": 0.7,
                    f"patch_icpp({3})%x_centroid": 0.9,  f"patch_icpp({3})%length_x": 0.2
                })

            if "x" in dimInfo[0]:
                dimParams[f"patch_icpp({patchID})%vel(1)"] = 0.0

            if "y" in dimInfo[0]:
                dimParams[f"patch_icpp({patchID})%vel(2)"] = 0.0

            if "z" in dimInfo[0]:
                dimParams[f"patch_icpp({patchID})%vel(3)"] = 0.0

        r.append((dimInfo, dimParams))

    return r


# pylint: disable=too-many-locals, too-many-statements
def generate_cases() -> typing.List[TestCase]:
    stack, cases = CaseGeneratorStack(), []

    def alter_bcs(dimInfo):
        for bc in [ -1, -2, -4, -5, -6, -7, -8, -9, -10, -11, -12, -3, -15, -16 ]:
            cases.append(create_case(stack, f"bc={bc}", get_bc_mods(bc, dimInfo)))

    def alter_weno():
        for weno_order in [3, 5]:
            stack.push(f"weno_order={weno_order}", {'weno_order': weno_order})

            for mapped_weno, mp_weno in [('F', 'F'), ('T', 'F'), ('F', 'T')]:
                trace = []
                if mapped_weno == 'T':
                    trace.append("mapped_weno")
                if mp_weno:
                    trace.append("mp_weno")

                if not (mp_weno == 'T' and weno_order != 5):
                    cases.append(create_case(stack, [
                        f"mapped_weno={mapped_weno}",
                        f"mp_weno={mp_weno}"
                    ], {
                        'mapped_weno': mapped_weno,
                        'mp_weno':     mp_weno
                    }))

            stack.pop()

    def alter_riemann_solvers(num_fluids):
        for riemann_solver in [1, 2]:
            stack.push(f"riemann_solver={riemann_solver}", {'riemann_solver': riemann_solver})

            cases.append(create_case(stack, "mixture_err",   {'mixture_err': 'T'}))
            cases.append(create_case(stack, "avg_state=1",   {'avg_state':   1}))
            cases.append(create_case(stack, "wave_speeds=2", {'wave_speeds': 2}))

            if riemann_solver == 2:
                cases.append(create_case(stack, "model_eqns=3", {'model_eqns': 3}))

            if num_fluids == 2:
                if riemann_solver == 2:
                    cases.append(create_case(stack, 'alt_soundspeed', {'alt_soundspeed': 'T'}))

                cases.append(create_case(stack, 'mpp_lim', {'mpp_lim': 'T'}))

            stack.pop()

    def alter_num_fluids(dimInfo):
        for num_fluids in [1, 2]:
            stack.push(f"{num_fluids} Fluid(s)", {"num_fluids": num_fluids})

            if num_fluids == 2:
                stack.push("", {
                    'fluid_pp(2)%gamma':          2.5,    'fluid_pp(2)%pi_inf':         0.0,  'patch_icpp(1)%alpha_rho(1)': 0.81,
                    'patch_icpp(1)%alpha(1)':     0.9,    'patch_icpp(1)%alpha_rho(2)': 0.19, 'patch_icpp(1)%alpha(2)':     0.1,
                    'patch_icpp(2)%alpha_rho(1)': 0.25,   'patch_icpp(2)%alpha(1)':     0.5,  'patch_icpp(2)%alpha_rho(2)': 0.25,
                    'patch_icpp(2)%alpha(2)':     0.5,    'patch_icpp(3)%alpha_rho(1)': 0.08, 'patch_icpp(3)%alpha(1)':     0.2,
                    'patch_icpp(3)%alpha_rho(2)': 0.0225, 'patch_icpp(3)%alpha(2)':     0.8
                })

            alter_riemann_solvers(num_fluids)
            alter_ib(dimInfo)

            if num_fluids == 1:
                stack.push("Viscous", {
                    'fluid_pp(1)%Re(1)' : 0.0001, 'dt' : 1e-11, 'patch_icpp(1)%vel(1)': 1.0})

                alter_ib(dimInfo, six_eqn_model=True)

                cases.append(create_case(stack, "",             {'weno_Re_flux': 'F'}))
                cases.append(create_case(stack, "weno_Re_flux", {'weno_Re_flux': 'T'}))

                for weno_Re_flux in ['T']:
                    stack.push("weno_Re_flux" if weno_Re_flux == 'T' else '', {'weno_Re_flux' : 'T'})
                    cases.append(create_case(stack, "weno_avg", {'weno_avg': 'T'}))
                    stack.pop()

                stack.pop()

            if num_fluids == 2:
                stack.push("Viscous", {
                    'fluid_pp(1)%Re(1)' : 0.001, 'fluid_pp(1)%Re(2)' : 0.001,
                    'fluid_pp(2)%Re(1)' : 0.001, 'fluid_pp(2)%Re(2)' : 0.001, 'dt' : 1e-11,
                    'patch_icpp(1)%vel(1)': 1.0}) 

                alter_ib(dimInfo, six_eqn_model=True)

                cases.append(create_case(stack, "",             {'weno_Re_flux': 'F'}))
                cases.append(create_case(stack, "weno_Re_flux", {'weno_Re_flux': 'T'}))
                for weno_Re_flux in ['T']:
                    stack.push("weno_Re_flux" if weno_Re_flux == 'T' else '', {'weno_Re_flux' : 'T'})
                    cases.append(create_case(stack, "weno_avg", {'weno_avg': 'T'}))
                    stack.pop()

                stack.pop()
                stack.pop()

            stack.pop()

    def alter_2d():
        stack.push("Axisymmetric", {
            'num_fluids' : 2, 'bc_y%beg': -2, 'cyl_coord': 'T',
            'fluid_pp(2)%gamma':          2.5,    'fluid_pp(2)%pi_inf':         0.0,  'patch_icpp(1)%alpha_rho(1)': 0.81,
            'patch_icpp(1)%alpha(1)':     0.9,    'patch_icpp(1)%alpha_rho(2)': 0.19, 'patch_icpp(1)%alpha(2)':     0.1,
            'patch_icpp(2)%alpha_rho(1)': 0.25,   'patch_icpp(2)%alpha(1)':     0.5,  'patch_icpp(2)%alpha_rho(2)': 0.25,
            'patch_icpp(2)%alpha(2)':     0.5,    'patch_icpp(3)%alpha_rho(1)': 0.08, 'patch_icpp(3)%alpha(1)':     0.2,
            'patch_icpp(3)%alpha_rho(2)': 0.0225, 'patch_icpp(3)%alpha(2)':     0.8,  'patch_icpp(1)%vel(1)': 0.0
        })

        cases.append(create_case(stack, "model_eqns=2", {'model_eqns': 2}))
        cases.append(create_case(stack, "model_eqns=3", {'model_eqns': 3}))

        stack.push("Viscous", {
            'fluid_pp(1)%Re(1)' : 0.0001, 'fluid_pp(1)%Re(2)' : 0.0001,
            'fluid_pp(2)%Re(1)' : 0.0001, 'fluid_pp(2)%Re(2)' : 0.0001, 'dt' : 1e-11}) 

        cases.append(create_case(stack, "",             {'weno_Re_flux': 'F'}))
        cases.append(create_case(stack, "weno_Re_flux", {'weno_Re_flux': 'T'}))
        for weno_Re_flux in ['T']:
            stack.push("weno_Re_flux" if weno_Re_flux == 'T' else '', {'weno_Re_flux' : 'T'})
            cases.append(create_case(stack, "weno_avg", {'weno_avg': 'T'}))
            stack.pop()

        stack.pop()
        stack.pop()

    def alter_3d():
        stack.push("Cylindrical", {
            'bc_y%beg': -14, 'bc_z%beg': -1, 'bc_z%end': -1, 'cyl_coord': 'T', 'x_domain%beg': 0.E+00, 
            'x_domain%end': 5.E+00, 'y_domain%beg': 0.E+00, 'y_domain%end': 1.E+00, 'z_domain%beg': 0.E+00,
            'z_domain%end' : 2.0*3.141592653589793E+00, 'm': 29, 'n': 29, 'p': 29,
            'patch_icpp(1)%geometry': 10, 'patch_icpp(1)%x_centroid' : 0.5, 'patch_icpp(1)%y_centroid' : 0.E+00,
            'patch_icpp(1)%z_centroid' : 0.E+00, 'patch_icpp(1)%radius' : 1.0, 'patch_icpp(1)%length_x' : 1.0,
            'patch_icpp(1)%length_y' : -1E+6, 'patch_icpp(1)%length_z' : -1E+6, 
            'patch_icpp(2)%geometry': 10, 'patch_icpp(2)%x_centroid' : 2.5, 'patch_icpp(2)%y_centroid' : 0.E+00,
            'patch_icpp(2)%z_centroid' : 0.E+00, 'patch_icpp(2)%radius' : 1.0, 'patch_icpp(2)%length_x' : 3.0,
            'patch_icpp(2)%length_y' : -1E+6, 'patch_icpp(2)%length_z' : -1E+6, 
            'patch_icpp(3)%geometry': 10, 'patch_icpp(3)%x_centroid' : 4.5, 'patch_icpp(3)%y_centroid' : 0.E+00,
            'patch_icpp(3)%z_centroid' : 0.E+00, 'patch_icpp(3)%radius' : 1.0, 'patch_icpp(3)%length_x' : 1.0,
            'patch_icpp(3)%length_y' : -1E+6, 'patch_icpp(3)%length_z' : -1E+6, 'patch_icpp(1)%vel(1)' :0.0,
             'num_fluids' : 2,
            'fluid_pp(2)%gamma':          2.5,    'fluid_pp(2)%pi_inf':         0.0,  'patch_icpp(1)%alpha_rho(1)': 0.81,
            'patch_icpp(1)%alpha(1)':     0.9,    'patch_icpp(1)%alpha_rho(2)': 0.19, 'patch_icpp(1)%alpha(2)':     0.1,
            'patch_icpp(2)%alpha_rho(1)': 0.25,   'patch_icpp(2)%alpha(1)':     0.5,  'patch_icpp(2)%alpha_rho(2)': 0.25,
            'patch_icpp(2)%alpha(2)':     0.5,    'patch_icpp(3)%alpha_rho(1)': 0.08, 'patch_icpp(3)%alpha(1)':     0.2,
            'patch_icpp(3)%alpha_rho(2)': 0.0225, 'patch_icpp(3)%alpha(2)':     0.8
        })

        cases.append(create_case(stack, "model_eqns=2", {'model_eqns': 2}))

        stack.push("Viscous", {
            'fluid_pp(1)%Re(1)' : 0.0001, 'fluid_pp(1)%Re(2)' : 0.0001,
            'fluid_pp(2)%Re(1)' : 0.0001, 'fluid_pp(2)%Re(2)' : 0.0001, 'dt' : 1e-11
            })

        cases.append(create_case(stack, "",             {'weno_Re_flux': 'F'}))
        cases.append(create_case(stack, "weno_Re_flux", {'weno_Re_flux': 'T'}))
        for weno_Re_flux in ['T']:
            stack.push("weno_Re_flux" if weno_Re_flux == 'T' else '', {'weno_Re_flux' : 'T'})
            cases.append(create_case(stack, "weno_avg", {'weno_avg': 'T'}))
            stack.pop()

        stack.pop()
        stack.pop()

    def alter_ppn(dimInfo):
        if len(dimInfo[0]) == 3:
            cases.append(create_case(stack, '2 MPI Ranks', {'m': 29, 'n': 29, 'p': 49}, ppn=2))
        else:
            cases.append(create_case(stack, '2 MPI Ranks', {}, ppn=2))

    def alter_ib(dimInfo, six_eqn_model=False):

        stack.push(f'IBM', {
            'ib': 'T', 'num_ibs': 1, 
            'patch_ib(1)%x_centroid': 0.5, 'patch_ib(1)%y_centroid': 0.5, 
            'patch_ib(1)%radius': 0.1, 'patch_icpp(1)%vel(1)': 0.001,
            'patch_icpp(2)%vel(1)': 0.001, 'patch_icpp(3)%vel(1)': 0.001,
        })

        if len(dimInfo[0]) == 3:
            cases.append(create_case(stack, f'', {
                'patch_ib(1)%z_centroid': 0.5,
                'patch_ib(1)%geometry': 8,
            }))
        elif len(dimInfo[0]) == 2:
            cases.append(create_case(stack, f'', {'patch_ib(1)%geometry': 2 }))
            if six_eqn_model:
                cases.append(create_case(stack, f'model_eqns=3', {'patch_ib(1)%geometry': 2, 'model_eqns': 3}))

        stack.pop()

    def alter_bubbles(dimInfo):
        if len(dimInfo[0]) > 0:
            stack.push("Bubbles", {"bubbles": 'T'})

            stack.push('', {
                'nb' : 3, 'fluid_pp(1)%gamma' : 0.16, 'fluid_pp(1)%pi_inf': 3515.0,
                'fluid_pp(2)%gamma': 2.5, 'fluid_pp(2)%pi_inf': 0.0, 'fluid_pp(1)%mul0' : 0.001002,
                'fluid_pp(1)%ss' : 0.07275,'fluid_pp(1)%pv' : 2338.8,'fluid_pp(1)%gamma_v' : 1.33,
                'fluid_pp(1)%M_v' : 18.02,'fluid_pp(1)%mu_v' : 8.816e-06,'fluid_pp(1)%k_v' : 0.019426,
                'fluid_pp(2)%gamma_v' : 1.4,'fluid_pp(2)%M_v' : 28.97,'fluid_pp(2)%mu_v' : 1.8e-05,
                'fluid_pp(2)%k_v' : 0.02556, 'patch_icpp(1)%alpha_rho(1)': 0.96, 'patch_icpp(1)%alpha(1)':
                4e-02, 'patch_icpp(2)%alpha_rho(1)': 0.96, 'patch_icpp(2)%alpha(1)': 4e-02,  'patch_icpp(3)%alpha_rho(1)': 0.96,
                'patch_icpp(3)%alpha(1)': 4e-02, 'patch_icpp(1)%pres': 1.0, 'patch_icpp(2)%pres': 1.0,
                'patch_icpp(3)%pres': 1.0
            })

            stack.push("Monopole", {"Monopole": 'T'})

            if len(dimInfo[0]) >= 2:
                stack.push("", {'Mono(1)%loc(2)': 0.5})

            if len(dimInfo[0]) >= 3:
                stack.push("", {'Mono(1)%loc(3)': 0.5, 'Mono(1)%support': 3})

            for polytropic in ['T', 'F']:
                stack.push("Polytropic" if polytropic == 'T' else '', {'polytropic' : polytropic})

                for bubble_model in [3, 2]:
                    stack.push(f"bubble_model={bubble_model}", {'bubble_model' : bubble_model})

                    if not (polytropic == 'F' and bubble_model == 3):
                        cases.append(create_case(stack, '', {}))

                    stack.pop()

                stack.pop()

            stack.push('', {'polytropic': 'T', 'bubble_model': 2})
            cases.append(create_case(stack, 'nb=1', {'nb': 1}))

            stack.push("QBMM", {'qbmm': 'T'})
            cases.append(create_case(stack, '', {}))


            stack.push("Non-polytropic", {'polytropic': 'F'})
            cases.append(create_case(stack, '', {}))

            stack.pop()

            stack.push('bubble_model=3', {'bubble_model': 3, 'polytropic': 'T'})
            cases.append(create_case(stack, '', {}))

            stack.push('Non-polytropic', { 'polytropic': 'F'})
            cases.append(create_case(stack, '', {}))

            for _ in range(7):
                stack.pop()

            if len(dimInfo[0]) >= 2:
                stack.pop()

            if len(dimInfo[0]) >= 3:
                stack.pop()

    def alter_hypoelasticity(dimInfo):
        # Hypoelasticity checks
        for num_fluids in [1,2]:
            stack.push(f"Hypoelasticity -> {num_fluids} Fluid(s)", {
                "hypoelasticity": 'T', "num_fluids": num_fluids,
                'riemann_solver':             1,
                'fluid_pp(1)%gamma':          0.3,    'fluid_pp(1)%pi_inf':         7.8E+05,
                'patch_icpp(1)%pres':         1.E+06, 'patch_icpp(1)%alpha_rho(1)': 1000.E+00,
                'patch_icpp(2)%pres':         1.E+05, 'patch_icpp(2)%alpha_rho(1)': 1000.E+00,
                'patch_icpp(3)%pres':         5.E+05, 'patch_icpp(3)%alpha_rho(1)': 1000.E+00,
                'patch_icpp(1)%tau_e(1)':     0.E+00, 'patch_icpp(2)%tau_e(1)':     0.E+00,
                'patch_icpp(3)%tau_e(1)':     0.E+00, 'fluid_pp(1)%G':              1.E+05,
            })

            if num_fluids == 2:
                stack.push("", {
                    'fluid_pp(2)%gamma':          0.3,    'fluid_pp(2)%pi_inf':      7.8E+05, 'patch_icpp(1)%alpha_rho(1)': 900.E+00,
                    'patch_icpp(1)%alpha(1)':     0.9,    'patch_icpp(1)%alpha_rho(2)':  100, 'patch_icpp(1)%alpha(2)':     0.1,
                    'patch_icpp(2)%alpha_rho(1)': 100,    'patch_icpp(2)%alpha(1)':     0.1,  'patch_icpp(2)%alpha_rho(2)': 900,
                    'patch_icpp(2)%alpha(2)':     0.9,    'patch_icpp(3)%alpha_rho(1)': 900,  'patch_icpp(3)%alpha(1)':     0.9,
                    'patch_icpp(3)%alpha_rho(2)': 100,    'patch_icpp(3)%alpha(2)':     0.1,
                    'fluid_pp(2)%G':              5.E+04
                })

            if len(dimInfo[0]) >= 2:
                stack.push("", {
                    'patch_icpp(1)%tau_e(2)':    0.E+00,  'patch_icpp(1)%tau_e(3)':     0.0E+00,
                    'patch_icpp(2)%tau_e(2)':    0.E+00,  'patch_icpp(2)%tau_e(3)':     0.0E+00,
                    'patch_icpp(3)%tau_e(2)':    0.E+00,  'patch_icpp(3)%tau_e(3)':     0.0E+00
                })

            if len(dimInfo[0]) == 3:
                stack.push("", {
                    'patch_icpp(1)%tau_e(4)': 0.E+00, 'patch_icpp(1)%tau_e(5)': 0.0E+00, 'patch_icpp(1)%tau_e(6)': 0.0E+00,
                    'patch_icpp(2)%tau_e(4)': 0.E+00, 'patch_icpp(2)%tau_e(5)': 0.0E+00, 'patch_icpp(2)%tau_e(6)': 0.0E+00,
                    'patch_icpp(3)%tau_e(4)': 0.E+00, 'patch_icpp(3)%tau_e(5)': 0.0E+00, 'patch_icpp(3)%tau_e(6)': 0.0E+00
                })

            cases.append(create_case(stack, '', {}))

            stack.pop()

            if num_fluids == 2:
                stack.pop()

            if len(dimInfo[0]) == 2:
                stack.pop()

            if len(dimInfo[0]) == 3:
                for _ in range(2):
                    stack.pop()

    def alter_phasechange(dimInfo):
        ndims = len(dimInfo[0])

        # Phase Change checks
        for relax_model in [5] + ([6] if ndims <= 2 else []):
            for num_fluids in ([2] if ndims == 1 or relax_model == 5 else []) + [3]:
                stack.push(f"Phase Change model {relax_model} -> {num_fluids} Fluid(s)", {
                    "relax": 'T',
                    "relax_model": relax_model,
                    'model_eqns': 3,
                    'palpha_eps': 1E-02,
                    'ptgalpha_eps': 1E-02,                    
                    "num_fluids": num_fluids,
                    'riemann_solver':           2,
                    'fluid_pp(1)%gamma':        0.7409,       'fluid_pp(1)%pi_inf':   1.7409E+09,
                    'fluid_pp(1)%cv':           1816,         'fluid_pp(1)%qv':   -1167000,
                    'fluid_pp(1)%qvp':          0.0,
                    'fluid_pp(2)%gamma':        2.3266,       'fluid_pp(2)%pi_inf':   0.0E+00,
                    'fluid_pp(2)%cv':           1040,         'fluid_pp(2)%qv':   2030000,
                    'fluid_pp(2)%qvp':          -23400,
                    'patch_icpp(1)%pres':       4.3755E+05,
                    'patch_icpp(1)%alpha(1)':   8.7149E-06,    'patch_icpp(1)%alpha_rho(1)': 9.6457E+02 * 8.7149E-06,
                    'patch_icpp(1)%alpha(2)':   1-8.7149E-06,  'patch_icpp(1)%alpha_rho(2)': 2.3132 * ( 1 - 8.7149E-06 ),
                    'patch_icpp(2)%pres':       9.6602E+04,
                    'patch_icpp(2)%alpha(1)':   3.6749E-05,   'patch_icpp(2)%alpha_rho(1)': 1.0957E+03 * 3.6749E-05,
                    'patch_icpp(2)%alpha(2)':   1-3.6749E-05, 'patch_icpp(2)%alpha_rho(2)': 0.5803 * ( 1 - 3.6749E-05 ),
                    'patch_icpp(3)%pres':       9.6602E+04,
                    'patch_icpp(3)%alpha(1)':   3.6749E-05,   'patch_icpp(3)%alpha_rho(1)': 1.0957E+03 * 3.6749E-05,
                    'patch_icpp(3)%alpha(2)':   1-3.6749E-05, 'patch_icpp(3)%alpha_rho(2)': 0.5803 * ( 1 - 3.6749E-05 )
                })

                if num_fluids == 3:
                    stack.push("", {
                        'fluid_pp(3)%gamma':        2.4870,        'fluid_pp(3)%pi_inf':   0.0E+00,
                        'fluid_pp(3)%cv':           717.5,         'fluid_pp(3)%qv':   0.0E+00,
                        'fluid_pp(3)%qvp':          0.0,
                        'patch_icpp(1)%alpha(2)':   2.5893E-02,                 'patch_icpp(1)%alpha_rho(2)':   2.3132 * 2.5893E-02,
                        'patch_icpp(2)%alpha(2)':   2.8728E-02,                 'patch_icpp(2)%alpha_rho(2)':   0.5803 * 2.8728E-02,
                        'patch_icpp(3)%alpha(2)':   2.8728E-02,                 'patch_icpp(3)%alpha_rho(2)':   0.5803 * 2.8728E-02,
                        'patch_icpp(1)%alpha(3)':   1-8.7149E-06-2.5893E-02,    'patch_icpp(1)%alpha_rho(3)':   3.5840 * ( 1-8.7149E-06-2.5893E-02 ),
                        'patch_icpp(2)%alpha(3)':   1-3.6749E-05-2.8728E-02,    'patch_icpp(2)%alpha_rho(3)':   0.8991 * ( 1-3.6749E-05-2.8728E-02 ),
                        'patch_icpp(3)%alpha(3)':   1-3.6749E-05-2.8728E-02,    'patch_icpp(3)%alpha_rho(3)':   0.8991 * ( 1-3.6749E-05-2.8728E-02 )
                    })

                if ndims == 1:
                    stack.push("", {
                        'patch_icpp(1)%vel(1)':   606.15, 'patch_icpp(2)%vel(1)': 10.0, 'patch_icpp(3)%vel(1)': 10.0
                    })
                elif ndims == 2:
                    stack.push("", {
                        'patch_icpp(1)%vel(1)':   0.0, 'patch_icpp(2)%vel(1)': 0.0, 'patch_icpp(3)%vel(1)': 0.0,
                        'patch_icpp(1)%vel(2)':   606.15, 'patch_icpp(2)%vel(2)': 10.0, 'patch_icpp(3)%vel(2)': 10.0
                    })
                elif ndims == 3:
                    stack.push("", {
                        'patch_icpp(1)%vel(1)':   0.0, 'patch_icpp(2)%vel(1)': 0.0, 'patch_icpp(3)%vel(1)': 0.0,
                        'patch_icpp(1)%vel(2)':   0.0, 'patch_icpp(2)%vel(2)': 0.0, 'patch_icpp(3)%vel(2)': 0.0,
                        'patch_icpp(1)%vel(3)':   606.15, 'patch_icpp(2)%vel(3)': 10.0, 'patch_icpp(3)%vel(3)': 10.0
                    })

                cases.append(create_case(stack, '', {}))

                stack.pop()
                stack.pop()

                if num_fluids == 3:
                    stack.pop()

    def alter_viscosity(dimInfo):
        # Viscosity & bubbles checks
        if len(dimInfo[0]) > 0:
            stack.push("Viscosity -> Bubbles",
                {"fluid_pp(1)%Re(1)": 50, "bubbles": 'T'})

            stack.push('', {
                'nb' : 1, 'fluid_pp(1)%gamma' : 0.16, 'fluid_pp(1)%pi_inf': 3515.0,
                'fluid_pp(2)%gamma': 2.5, 'fluid_pp(2)%pi_inf': 0.0, 'fluid_pp(1)%mul0' : 0.001002,
                'fluid_pp(1)%ss' : 0.07275,'fluid_pp(1)%pv' : 2338.8,'fluid_pp(1)%gamma_v' : 1.33,
                'fluid_pp(1)%M_v' : 18.02,'fluid_pp(1)%mu_v' : 8.816e-06,'fluid_pp(1)%k_v' : 0.019426,
                'fluid_pp(2)%gamma_v' : 1.4,'fluid_pp(2)%M_v' : 28.97,'fluid_pp(2)%mu_v' : 1.8e-05,
                'fluid_pp(2)%k_v' : 0.02556, 'patch_icpp(1)%alpha_rho(1)': 0.96, 'patch_icpp(1)%alpha(1)': 4e-02,
                'patch_icpp(2)%alpha_rho(1)': 0.96, 'patch_icpp(2)%alpha(1)': 4e-02,  'patch_icpp(3)%alpha_rho(1)': 0.96,
                'patch_icpp(3)%alpha(1)': 4e-02, 'patch_icpp(1)%pres': 1.0, 'patch_icpp(2)%pres': 1.0,
                'patch_icpp(3)%pres': 1.0
            })

            for polytropic in ['T', 'F']:
                stack.push("Polytropic" if polytropic == 'T' else '', {'polytropic' : polytropic})

                for bubble_model in [3, 2]:
                    stack.push(f"bubble_model={bubble_model}", {'bubble_model' : bubble_model})

                    if not (polytropic == 'F' and bubble_model == 3):
                        cases.append(create_case(stack, '', {}))

                    stack.pop()

                stack.pop()

            stack.push('', {'polytropic': 'T', 'bubble_model': 2})
            cases.append(create_case(stack, 'nb=1', {'nb': 1}))

            stack.push("QBMM", {'qbmm': 'T'})
            cases.append(create_case(stack, '', {}))

            stack.push('bubble_model=3', {'bubble_model': 3})
            cases.append(create_case(stack, '', {}))

            for _ in range(5):
                stack.pop()

    def foreach_dimension():
        for dimInfo, dimParams in get_dimensions():
            stack.push(f"{len(dimInfo[0])}D", dimParams)
            alter_bcs(dimInfo)
            alter_weno()
            alter_num_fluids(dimInfo)
            if len(dimInfo[0]) == 2:
                alter_2d()
            if len(dimInfo[0]) == 3:
                alter_3d()
            alter_ppn(dimInfo)
            stack.push('', {'dt': [1e-07, 1e-06, 1e-06][len(dimInfo[0])-1]})
            alter_bubbles(dimInfo)
            alter_hypoelasticity(dimInfo)
            alter_phasechange(dimInfo)
            alter_viscosity(dimInfo)
            stack.pop()
            stack.pop()

    foreach_dimension()

    # Sanity Check 1
    if stack.size() != 0:
        raise common.MFCException("generate_cases: stack isn't fully pop'ed")

    # Sanity Check 2
    uuids  = [ case.get_uuid() for case in cases ]
    l1, l2 = len(uuids), len(set(uuids))
    if l1 != l2:
        raise common.MFCException(f"generate_cases: uuids aren't unique ({l1} cases but {l2} unique uuids)")

    return cases
