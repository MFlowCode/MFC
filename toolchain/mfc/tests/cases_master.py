import common

from tests.case import Case
from tests.case import create_case
from tests.case import CaseGeneratorStack


def get_dimensions():
    r = []

    for dimInfo in [ (["x"],           {'m': 299, 'n': 0,  'p': 0},  {"geometry": 1}),
                     (["x", "y"],      {'m': 49,  'n': 39, 'p': 0},  {"geometry": 3}),
                     (["x", "y", "z"], {'m': 24,  'n': 24, 'p': 24}, {"geometry": 9}) ]:
        dimParams = {**dimInfo[1]}

        for dimCmp in dimInfo[0]:
            dimParams.update({ f"{dimCmp}_domain%beg": 0.E+00,
                               f"{dimCmp}_domain%end": 1.E+00 })

        for patchID in range(1, 3+1):
            dimParams[f"patch_icpp({patchID})%geometry"] = dimInfo[2].get("geometry")

            if "z" in dimInfo[0]:
                dimParams.update({
                    f"patch_icpp({1})%z_centroid":       0.05,
                    f"patch_icpp({1})%length_z":         0.1,
                    f"patch_icpp({2})%z_centroid":       0.45,
                    f"patch_icpp({2})%length_z":         0.7,
                    f"patch_icpp({3})%z_centroid":       0.9,
                    f"patch_icpp({3})%length_z":         0.2,
                    f"patch_icpp({patchID})%y_centroid": 0.5,
                    f"patch_icpp({patchID})%length_y":   1,
                    f"patch_icpp({patchID})%x_centroid": 0.5,
                    f"patch_icpp({patchID})%length_x":   1
                })

            elif "y" in dimInfo[0]:
                dimParams.update({
                    f"patch_icpp({1})%y_centroid": 0.05,
                    f"patch_icpp({1})%length_y": 0.1,
                    f"patch_icpp({2})%y_centroid": 0.45,
                    f"patch_icpp({2})%length_y": 0.7,
                    f"patch_icpp({3})%y_centroid": 0.9,
                    f"patch_icpp({3})%length_y": 0.2,
                    f"patch_icpp({patchID})%x_centroid": 0.5,
                    f"patch_icpp({patchID})%length_x": 1
                })
            else:
                dimParams.update({
                    f"patch_icpp({1})%x_centroid": 0.05,
                    f"patch_icpp({1})%length_x":   0.1,
                    f"patch_icpp({2})%x_centroid": 0.45,
                    f"patch_icpp({2})%length_x":   0.7,
                    f"patch_icpp({3})%x_centroid": 0.9,
                    f"patch_icpp({3})%length_x":   0.2
                })

            if "x" in dimInfo[0]:
                dimParams[f"patch_icpp({patchID})%vel(1)"] = 0.0

            if "y" in dimInfo[0]:
                dimParams[f"patch_icpp({patchID})%vel(2)"] = 0.0

            if "z" in dimInfo[0]:
                dimParams[f"patch_icpp({patchID})%vel(3)"] = 0.0

        r.append((dimInfo, dimParams))

    return r


def generate_cases() -> list:
    stack, cases = CaseGeneratorStack(), []

    for dimInfo, dimParams in get_dimensions():
        stack.push(f"{len(dimInfo[0])}D (m={dimInfo[1].get('m')},n={dimInfo[1].get('n')},p={dimInfo[1].get('p')})", dimParams)

        for bc in [ -1, -2, -4, -5, -6, -7, -8, -9, -10, -11, -12, -3 ]:
            params = {}
            for dimCmp in dimInfo[0]:
                params.update({f'bc_{dimCmp}%beg': bc, f'bc_{dimCmp}%end': bc})

            stack.push(f"bc={bc}", params)
            cases.append(create_case(stack, '', {}))

            if bc != -3: # Use bc = 3 henceforth
                stack.pop()

        for weno_order in [3, 5]:
            stack.push(f"weno_order={weno_order}", {'weno_order': weno_order})

            for mapped_weno, mp_weno in [('F', 'F'), ('T', 'F'), ('F', 'T')]:
                stack.push(f"(mapped_weno={mapped_weno},mp_weno={mp_weno})", {
                    'mapped_weno': mapped_weno,
                    'mp_weno':     mp_weno
                })

                if not (mp_weno == 'T' and weno_order != 5):
                    cases.append(create_case(stack, '', {}))

                stack.pop()

            stack.pop()

        for num_fluids in [1, 2]:
            stack.push(f"num_fluids={num_fluids}", {"num_fluids": num_fluids})

            if num_fluids == 2:
                stack.push("", {
                    'fluid_pp(2)%gamma':          2.5,    'fluid_pp(2)%pi_inf':         0.0,  'patch_icpp(1)%alpha_rho(1)': 0.81,
                    'patch_icpp(1)%alpha(1)':     0.9,    'patch_icpp(1)%alpha_rho(2)': 0.19, 'patch_icpp(1)%alpha(2)':     0.1,
                    'patch_icpp(2)%alpha_rho(1)': 0.25,   'patch_icpp(2)%alpha(1)':     0.5,  'patch_icpp(2)%alpha_rho(2)': 0.25,
                    'patch_icpp(2)%alpha(2)':     0.5,    'patch_icpp(3)%alpha_rho(1)': 0.08, 'patch_icpp(3)%alpha(1)':     0.2,
                    'patch_icpp(3)%alpha_rho(2)': 0.0225, 'patch_icpp(3)%alpha(2)':     0.8
                })

            for riemann_solver in [1, 2]:
                stack.push(f"riemann_solver={riemann_solver}", {'riemann_solver': riemann_solver})

                cases.append(create_case(stack, "mixture_err=T", {'mixture_err': 'T'}))
                cases.append(create_case(stack, "avg_state=1",   {'avg_state':   '1'}))
                cases.append(create_case(stack, "wave_speeds=2", {'wave_speeds': '2'}))
                if riemann_solver == 2:
                    cases.append(create_case(stack, "model_eqns=3", {'model_eqns': 3}))

                if num_fluids == 2:
                    if riemann_solver == 2:
                        cases.append(create_case(stack, 'alt_soundspeed=T', {'alt_soundspeed': 'T'}))

                    cases.append(create_case(stack, 'mpp_lim=T', {'mpp_lim': 'T'}))

                stack.pop()

            if num_fluids == 2:
                stack.pop()

            stack.pop()

        for num_fluids in [2]:
            stack.push(f"num_fluids={num_fluids}", {"num_fluids": num_fluids})


            stack.push("", {
                'fluid_pp(2)%gamma':          2.5,    'fluid_pp(2)%pi_inf':         0.0,  'patch_icpp(1)%alpha_rho(1)': 0.81,
                'patch_icpp(1)%alpha(1)':     0.9,    'patch_icpp(1)%alpha_rho(2)': 0.19, 'patch_icpp(1)%alpha(2)':     0.1,
                'patch_icpp(2)%alpha_rho(1)': 0.25,   'patch_icpp(2)%alpha(1)':     0.5,  'patch_icpp(2)%alpha_rho(2)': 0.25,
                'patch_icpp(2)%alpha(2)':     0.5,    'patch_icpp(3)%alpha_rho(1)': 0.08, 'patch_icpp(3)%alpha(1)':     0.2,
                'patch_icpp(3)%alpha_rho(2)': 0.0225, 'patch_icpp(3)%alpha(2)':     0.8
            })

            stack.push(f"Viscous", {'fluid_pp(1)%Re(1)' : 0.001, 'fluid_pp(1)%Re(2)' : 0.001,
                'fluid_pp(2)%Re(1)' : 0.001, 'fluid_pp(2)%Re(2)' : 0.001, 'dt' : 1e-11}) 

            cases.append(create_case(stack, "weno_Re_flux=F", {'weno_Re_flux': 'F'}))             
            cases.append(create_case(stack, "weno_Re_flux=T", {'weno_Re_flux': 'T'}))

            stack.pop() 

            stack.pop()

            stack.pop()

        if len(dimInfo[0]) == 2:
            stack.push(f"Axisymmetric", {'bc_y%beg': -2, 'cyl_coord': 'T'})

            stack.push("", { 'num_fluids' : 2,
                'fluid_pp(2)%gamma':          2.5,    'fluid_pp(2)%pi_inf':         0.0,  'patch_icpp(1)%alpha_rho(1)': 0.81,
                'patch_icpp(1)%alpha(1)':     0.9,    'patch_icpp(1)%alpha_rho(2)': 0.19, 'patch_icpp(1)%alpha(2)':     0.1,
                'patch_icpp(2)%alpha_rho(1)': 0.25,   'patch_icpp(2)%alpha(1)':     0.5,  'patch_icpp(2)%alpha_rho(2)': 0.25,
                'patch_icpp(2)%alpha(2)':     0.5,    'patch_icpp(3)%alpha_rho(1)': 0.08, 'patch_icpp(3)%alpha(1)':     0.2,
                'patch_icpp(3)%alpha_rho(2)': 0.0225, 'patch_icpp(3)%alpha(2)':     0.8
            })

            cases.append(create_case(stack, "model_eqns=2", {'model_eqns': 2}))
            cases.append(create_case(stack, "model_eqns=3", {'model_eqns': 3}))

            stack.push(f"Viscous", {'fluid_pp(1)%Re(1)' : 0.0001, 'fluid_pp(1)%Re(2)' : 0.0001,
                'fluid_pp(2)%Re(1)' : 0.0001, 'fluid_pp(2)%Re(2)' : 0.0001, 'dt' : 1e-11}) 

            cases.append(create_case(stack, "weno_Re_flux=F", {'weno_Re_flux': 'F'}))             
            cases.append(create_case(stack, "weno_Re_flux=T", {'weno_Re_flux': 'T'}))

            stack.pop()

            stack.pop()

            stack.pop()  

        if(len(dimInfo[0]) == 3):
            stack.push(f"Cylindrical", {'bc_y%beg': -13, 'bc_z%beg': -1, 'bc_z%end': -1, 'cyl_coord': 'T', 'x_domain%beg': 0.E+00, 
                'x_domain%end': 5.E+00, 'y_domain%beg': 0.E+00, 'y_domain%end': 1.E+00, 'z_domain%beg': 0.E+00, 'z_domain%end' : 2.0*3.141592653589793E+00,
                'm': 29, 'n': 29, 'p': 29})


            stack.push("", { 'patch_icpp(1)%geometry': 10, 'patch_icpp(1)%x_centroid' : 0.5, 'patch_icpp(1)%y_centroid' : 0.E+00,
                'patch_icpp(1)%z_centroid' : 0.E+00, 'patch_icpp(1)%radius' : 1.0, 'patch_icpp(1)%length_x' : 1.0,
                'patch_icpp(1)%length_y' : -1E+6, 'patch_icpp(1)%length_z' : -1E+6, 
                'patch_icpp(2)%geometry': 10, 'patch_icpp(2)%x_centroid' : 2.5, 'patch_icpp(2)%y_centroid' : 0.E+00,
                'patch_icpp(2)%z_centroid' : 0.E+00, 'patch_icpp(2)%radius' : 1.0, 'patch_icpp(2)%length_x' : 3.0,
                'patch_icpp(2)%length_y' : -1E+6, 'patch_icpp(2)%length_z' : -1E+6, 
                'patch_icpp(3)%geometry': 10, 'patch_icpp(3)%x_centroid' : 4.5, 'patch_icpp(3)%y_centroid' : 0.E+00,
                'patch_icpp(3)%z_centroid' : 0.E+00, 'patch_icpp(3)%radius' : 1.0, 'patch_icpp(3)%length_x' : 1.0,
                'patch_icpp(3)%length_y' : -1E+6, 'patch_icpp(3)%length_z' : -1E+6})

            stack.push("", { 'num_fluids' : 2,
                'fluid_pp(2)%gamma':          2.5,    'fluid_pp(2)%pi_inf':         0.0,  'patch_icpp(1)%alpha_rho(1)': 0.81,
                'patch_icpp(1)%alpha(1)':     0.9,    'patch_icpp(1)%alpha_rho(2)': 0.19, 'patch_icpp(1)%alpha(2)':     0.1,
                'patch_icpp(2)%alpha_rho(1)': 0.25,   'patch_icpp(2)%alpha(1)':     0.5,  'patch_icpp(2)%alpha_rho(2)': 0.25,
                'patch_icpp(2)%alpha(2)':     0.5,    'patch_icpp(3)%alpha_rho(1)': 0.08, 'patch_icpp(3)%alpha(1)':     0.2,
                'patch_icpp(3)%alpha_rho(2)': 0.0225, 'patch_icpp(3)%alpha(2)':     0.8
            })

            cases.append(create_case(stack, "model_eqns=2", {'model_eqns': 2}))

            stack.push(f"Viscous", {'fluid_pp(1)%Re(1)' : 0.0001, 'fluid_pp(1)%Re(2)' : 0.0001,
                'fluid_pp(2)%Re(1)' : 0.0001, 'fluid_pp(2)%Re(2)' : 0.0001, 'dt' : 1e-11}) 

            cases.append(create_case(stack, "weno_Re_flux=F", {'weno_Re_flux': 'F'}))             
            cases.append(create_case(stack, "weno_Re_flux=T", {'weno_Re_flux': 'T'}))

            stack.pop()

            stack.pop()

            stack.pop()

            stack.pop() 




        if len(dimInfo[0]) == 3:
            cases.append(create_case(stack, f'ppn=2,m=29,n=29,p=49', {'m': 29, 'n': 29, 'p': 49}, ppn=2))
        else:
            cases.append(create_case(stack, f'ppn=2', {}, ppn=2))

        stack.push('', {'dt': [1e-07, 1e-06, 1e-06][len(dimInfo[0])-1]})

        if len(dimInfo[0]) > 0:
            stack.push(f"bubbles={'T'}", {"bubbles": 'T'})

            stack.push(f'', {
                'nb' : 3, 'fluid_pp(1)%gamma' : 0.16, 'fluid_pp(1)%pi_inf': 3515.0,
                'fluid_pp(2)%gamma': 2.5, 'fluid_pp(2)%pi_inf': 0.0, 'fluid_pp(1)%mul0' : 0.001002,
                'fluid_pp(1)%ss' : 0.07275,'fluid_pp(1)%pv' : 2338.8,'fluid_pp(1)%gamma_v' : 1.33,
                'fluid_pp(1)%M_v' : 18.02,'fluid_pp(1)%mu_v' : 8.816e-06,'fluid_pp(1)%k_v' : 0.019426,
                'fluid_pp(2)%gamma_v' : 1.4,'fluid_pp(2)%M_v' : 28.97,'fluid_pp(2)%mu_v' : 1.8e-05,
                'fluid_pp(2)%k_v' : 0.02556, 'patch_icpp(1)%alpha_rho(1)': 0.999999999999, 'patch_icpp(1)%alpha(1)':
                1e-12, 'patch_icpp(2)%alpha_rho(1)': 0.96, 'patch_icpp(2)%alpha(1)': 4e-02,  'patch_icpp(3)%alpha_rho(1)': 0.999999999999,
                'patch_icpp(3)%alpha(1)': 1e-12, 'patch_icpp(1)%pres': 1.0, 'patch_icpp(2)%pres': 1.0,
                'patch_icpp(3)%pres': 1.0
            })

            stack.push(f"Monopole={'T'}", {"Monopole": 'T'})

            if len(dimInfo[0]) >= 2:
                stack.push("", {'Mono(1)%loc(2)': 0.5})

            if len(dimInfo[0]) >= 3:
                stack.push("", {'Mono(1)%loc(3)': 0.5, 'Mono(1)%support': 3})

            for polytropic in ['T', 'F']:
                stack.push(f"polytropic={polytropic}", {'polytropic' : polytropic})

                for bubble_model in [3, 2]:
                    stack.push(f"bubble_model={bubble_model}", {'bubble_model' : bubble_model})

                    if not (polytropic == 'F' and bubble_model == 3):
                        cases.append(create_case(stack, '', {}))

                    stack.pop()

                stack.pop()

            stack.push('', {'polytropic': 'T', 'bubble_model': 2})
            cases.append(create_case(stack, 'nb=1', {'nb': 1}))

            stack.push(f"qbmm={'T'}", {'qbmm': 'T'})
            cases.append(create_case(stack, '', {}))

            stack.push('bubble_model=3', {'bubble_model': 3})
            cases.append(create_case(stack, '', {}))

            for i in range(6):
                stack.pop()

            if len(dimInfo[0]) >= 2:
                stack.pop()

            if len(dimInfo[0]) >= 3:
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


def generate_filtered_cases(args: dict):
    cases: list = generate_cases()

    if len(args["only"]) > 0:
        for i, case in enumerate(cases[:]):
            case: Case

            doKeep = False
            for o in args["only"]:
                if str(o) == case.get_uuid():
                    doKeep = True
                    break

            if not doKeep:
                cases.remove(case)

    return cases
