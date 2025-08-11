# pylint: disable=too-many-lines
import os, typing, itertools

from mfc   import common
from .case import Nt, define_case_d, define_case_f, CaseGeneratorStack, TestCaseBuilder

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
def list_cases() -> typing.List[TestCaseBuilder]:
    stack, cases = CaseGeneratorStack(), []

    def alter_bcs(dimInfo):
        for bc in [ -1, -2, -4, -5, -6, -7, -8, -9, -10, -11, -12, -3, -15, -16, -17]:
            cases.append(define_case_d(stack, f"bc={bc}", get_bc_mods(bc, dimInfo)))

    def alter_grcbc(dimInfo):
        if len(dimInfo[0]) == 1:
            stack.push('', {'patch_icpp(1)%vel(1)':1.0, 'patch_icpp(2)%vel(1)':1.0, 'patch_icpp(3)%vel(1)':1.0,
                            'bc_x%beg':-7, 'bc_x%end':-8, 'bc_x%grcbc_in':'T', 'bc_x%grcbc_out':'T', 'bc_x%grcbc_vel_out':'T',
                            'bc_x%vel_in(1)':1.0, 'bc_x%vel_in(2)':0.0, 'bc_x%vel_in(3)':0.0, 'bc_x%vel_out(1)':1.0, 'bc_x%vel_out(2)':0.0, 'bc_x%vel_out(3)':0.0,
                            'bc_x%pres_in':1.0, 'bc_x%pres_out':1.0, 'bc_x%alpha_in(1)':1.0, 'bc_x%alpha_rho_in(1)':1.0})
            cases.append(define_case_d(stack, [f"grcbc x"],{}))
            stack.pop()
        elif len(dimInfo[0]) == 2:
            stack.push('', {'patch_icpp(1)%vel(1)':1.0, 'patch_icpp(2)%vel(1)':1.0, 'patch_icpp(3)%vel(1)':1.0,
                            'bc_x%beg':-7, 'bc_x%end':-8, 'bc_x%grcbc_in':'T', 'bc_x%grcbc_out':'T', 'bc_x%grcbc_vel_out':'T',
                            'bc_x%vel_in(1)':1.0, 'bc_x%vel_in(2)':0.0, 'bc_x%vel_in(3)':0.0, 'bc_x%vel_out(1)':1.0, 'bc_x%vel_out(2)':0.0, 'bc_x%vel_out(3)':0.0,
                            'bc_x%pres_in':1.0, 'bc_x%pres_out':1.0, 'bc_x%alpha_in(1)':1.0, 'bc_x%alpha_rho_in(1)':1.0})
            cases.append(define_case_d(stack, [f"grcbc x"],{}))
            stack.pop()

            stack.push('', {'patch_icpp(1)%vel(2)':1.0, 'patch_icpp(2)%vel(2)':1.0, 'patch_icpp(3)%vel(2)':1.0,
                            'bc_y%beg':-7, 'bc_y%end':-8, 'bc_y%grcbc_in':'T', 'bc_y%grcbc_out':'T', 'bc_y%grcbc_vel_out':'T',
                            'bc_y%vel_in(1)':0.0, 'bc_y%vel_in(2)':1.0, 'bc_y%vel_in(3)':0.0, 'bc_y%vel_out(1)':0.0, 'bc_y%vel_out(2)':1.0, 'bc_y%vel_out(3)':0.0,
                            'bc_y%pres_in':1.0, 'bc_y%pres_out':1.0, 'bc_y%alpha_in(1)':1.0, 'bc_y%alpha_rho_in(1)':1.0})
            cases.append(define_case_d(stack, [f"grcbc y"],{}))
            stack.pop()
        elif len(dimInfo[0]) == 3:
            stack.push('', {'patch_icpp(1)%vel(1)':1.0, 'patch_icpp(2)%vel(1)':1.0, 'patch_icpp(3)%vel(1)':1.0,
                            'bc_x%beg':-7, 'bc_x%end':-8, 'bc_x%grcbc_in':'T', 'bc_x%grcbc_out':'T', 'bc_x%grcbc_vel_out':'T',
                            'bc_x%vel_in(1)':1.0, 'bc_x%vel_in(2)':0.0, 'bc_x%vel_in(3)':0.0, 'bc_x%vel_out(1)':1.0, 'bc_x%vel_out(2)':0.0, 'bc_x%vel_out(3)':0.0,
                            'bc_x%pres_in':1.0, 'bc_x%pres_out':1.0, 'bc_x%alpha_in(1)':1.0, 'bc_x%alpha_rho_in(1)':1.0})
            cases.append(define_case_d(stack, [f"grcbc x"],{}))
            stack.pop()

            stack.push('', {'patch_icpp(1)%vel(2)':1.0, 'patch_icpp(2)%vel(2)':1.0, 'patch_icpp(3)%vel(2)':1.0,
                            'bc_y%beg':-7, 'bc_y%end':-8, 'bc_y%grcbc_in':'T', 'bc_y%grcbc_out':'T', 'bc_y%grcbc_vel_out':'T',
                            'bc_y%vel_in(1)':0.0, 'bc_y%vel_in(2)':1.0, 'bc_y%vel_in(3)':0.0, 'bc_y%vel_out(1)':0.0, 'bc_y%vel_out(2)':1.0, 'bc_y%vel_out(3)':0.0,
                            'bc_y%pres_in':1.0, 'bc_y%pres_out':1.0, 'bc_y%alpha_in(1)':1.0, 'bc_y%alpha_rho_in(1)':1.0})
            cases.append(define_case_d(stack, [f"grcbc y"],{}))
            stack.pop()

            stack.push('', {'patch_icpp(1)%vel(3)':1.0, 'patch_icpp(2)%vel(3)':1.0, 'patch_icpp(3)%vel(3)':1.0,
                            'bc_z%beg':-7, 'bc_z%end':-8, 'bc_z%grcbc_in':'T', 'bc_z%grcbc_out':'T', 'bc_z%grcbc_vel_out':'T',
                            'bc_z%vel_in(1)':0.0, 'bc_z%vel_in(2)':0.0, 'bc_z%vel_in(3)':1.0, 'bc_z%vel_out(1)':0.0, 'bc_z%vel_out(2)':0.0, 'bc_z%vel_out(3)':1.0,
                            'bc_z%pres_in':1.0, 'bc_z%pres_out':1.0, 'bc_z%alpha_in(1)':1.0, 'bc_z%alpha_rho_in(1)':1.0})
            cases.append(define_case_d(stack, [f"grcbc z"],{}))
            stack.pop()

    def alter_capillary():
        stack.push('', {'patch_icpp(1)%cf_val':1, 'patch_icpp(2)%cf_val':0, 'patch_icpp(3)%cf_val':1,
                        'sigma':1, 'model_eqns':3, 'surface_tension': 'T'})
        cases.append(define_case_d(stack, [f"capillary=T","model_eqns=3"],{}))
        stack.pop()

    def alter_weno(dimInfo):
        for weno_order in [3, 5, 7]:
            stack.push(f"weno_order={weno_order}", {'weno_order': weno_order})
            for mapped_weno, wenoz, teno, mp_weno in itertools.product('FT', repeat=4):

                if sum(var == 'T' for var in [mapped_weno, wenoz, teno, mp_weno]) > 1:
                    continue
                if mp_weno == 'T' and weno_order != 5:
                    continue
                if teno == 'T' and weno_order == 3:
                    continue

                trace = [f"{var}={val}" for var, val in zip(["mapped_weno", "wenoz", "teno", "mp_weno"], [mapped_weno, wenoz, teno, mp_weno]) if val == 'T']
                data = {var: 'T' for var, val in zip(["mapped_weno", "wenoz", "teno", "mp_weno"], [mapped_weno, wenoz, teno, mp_weno]) if val == 'T'}

                if "teno" in data:
                    data["teno_CT"] = 1e-6
                if "wenoz" in data and weno_order == 7:
                    data["wenoz_q"] = 3.0

                if weno_order == 7:
                    data = {**data, 'weno_eps': 1e-6} # increase damping for stability

                    if "z" in dimInfo[0]:
                        data = {**data, 'm': 35, 'n': 35, 'p': 35}

                cases.append(define_case_d(stack, trace, data))

            stack.pop()

    def alter_igr():
        stack.push('IGR',{'igr': 'T',  'alf_factor': 10, 'num_igr_iters': 10,
                       'elliptic_smoothing': 'T', 'elliptic_smoothing_iters': 10,
                       'num_igr_warm_start_iters': 10})

        for order in [3, 5]:
            stack.push(f"igr_order={order}", {'igr_order': order})

            cases.append(define_case_d(stack, 'Jacobi', {'igr_iter_solver': 1}))
            if order == 5:
                cases.append(define_case_d(stack, 'Gauss Seidel', {'igr_iter_solver': 2}))

            stack.pop()

        stack.pop()

    def alter_muscl():
        for muscl_order in [1, 2]:
            stack.push(f"muscl_order={muscl_order}", {'muscl_order': muscl_order, 'recon_type':2, 'weno_order':0})

            if muscl_order == 1:
                for int_comp in ["T", "F"]:
                    cases.append(define_case_d(stack, f"int_comp={int_comp}", {'int_comp': int_comp}))
            elif muscl_order == 2:
                for int_comp in ["T", "F"]:
                    stack.push(f"int_comp={int_comp}", {'int_comp': int_comp})
                    cases.append(define_case_d(stack, f"muscl_lim=1", {'muscl_lim': 1}))
                    stack.pop()
                for muscl_lim in [2,3,4,5]:
                    cases.append(define_case_d(stack, f"muscl_lim={muscl_lim}", {'muscl_lim': muscl_lim}))
            stack.pop()

    def alter_riemann_solvers(num_fluids):
        for riemann_solver in [1, 2]:
            stack.push(f"riemann_solver={riemann_solver}", {'riemann_solver': riemann_solver})

            cases.append(define_case_d(stack, "mixture_err",   {'mixture_err': 'T'}))
            cases.append(define_case_d(stack, "avg_state=1",   {'avg_state':   1}))
            cases.append(define_case_d(stack, "wave_speeds=2", {'wave_speeds': 2}))

            if riemann_solver == 2:
                cases.append(define_case_d(stack, "model_eqns=3", {'model_eqns': 3}))

            if num_fluids == 2:
                if riemann_solver == 2:
                    cases.append(define_case_d(stack, 'alt_soundspeed', {'alt_soundspeed': 'T'}))

                cases.append(define_case_d(stack, 'mpp_lim', {'mpp_lim': 'T'}))

            stack.pop()

    def alter_low_Mach_correction():
        stack.push('', {'fluid_pp(1)%gamma' : 0.16, 'fluid_pp(1)%pi_inf': 3515.0, 'dt': 1e-7})

        stack.push(f"riemann_solver=1",{'riemann_solver': 1})
        cases.append(define_case_d(stack, 'low_Mach=1', {'low_Mach': 1}))
        stack.pop()
        stack.push(f"riemann_solver=2",{'riemann_solver': 2})
        cases.append(define_case_d(stack, 'low_Mach=1', {'low_Mach': 1}))
        cases.append(define_case_d(stack, 'low_Mach=2', {'low_Mach': 2}))
        stack.pop()

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

                if len(dimInfo[0]) > 1:
                    alter_capillary()

            alter_riemann_solvers(num_fluids)
            alter_low_Mach_correction()
            alter_ib(dimInfo)
            if len(dimInfo[0]) > 1:
                alter_igr()


            if num_fluids == 1:

                stack.push("Viscous", {
                    'fluid_pp(1)%Re(1)' : 0.0001, 'dt' : 1e-11, 'patch_icpp(1)%vel(1)': 1.0,
                    'viscous': 'T'})

                alter_ib(dimInfo, six_eqn_model=True)

                if len(dimInfo[0]) > 1:
                    alter_igr()

                cases.append(define_case_d(stack, "",             {'weno_Re_flux': 'F'}))
                cases.append(define_case_d(stack, "weno_Re_flux", {'weno_Re_flux': 'T'}))

                for weno_Re_flux in ['T']:
                    stack.push("weno_Re_flux" if weno_Re_flux == 'T' else '', {'weno_Re_flux' : 'T'})
                    cases.append(define_case_d(stack, "weno_avg", {'weno_avg': 'T'}))
                    stack.pop()

                stack.pop()

            if num_fluids == 2:
                stack.push("Viscous", {
                    'fluid_pp(1)%Re(1)' : 0.001, 'fluid_pp(1)%Re(2)' : 0.001,
                    'fluid_pp(2)%Re(1)' : 0.001, 'fluid_pp(2)%Re(2)' : 0.001, 'dt' : 1e-11,
                    'patch_icpp(1)%vel(1)': 1.0, 'viscous': 'T'})

                alter_ib(dimInfo, six_eqn_model=True)

                if len(dimInfo[0]) > 1:
                    alter_igr()

                cases.append(define_case_d(stack, "",             {'weno_Re_flux': 'F'}))
                cases.append(define_case_d(stack, "weno_Re_flux", {'weno_Re_flux': 'T'}))
                for weno_Re_flux in ['T']:
                    stack.push("weno_Re_flux" if weno_Re_flux == 'T' else '', {'weno_Re_flux' : 'T'})
                    cases.append(define_case_d(stack, "weno_avg", {'weno_avg': 'T'}))
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

        cases.append(define_case_d(stack, "model_eqns=2", {'model_eqns': 2}))
        cases.append(define_case_d(stack, "model_eqns=3", {'model_eqns': 3}))
        cases.append(define_case_d(stack, "HLL", {'riemann_solver': 1}))

        stack.push("Viscous", {
            'fluid_pp(1)%Re(1)' : 0.0001, 'fluid_pp(1)%Re(2)' : 0.0001,
            'fluid_pp(2)%Re(1)' : 0.0001, 'fluid_pp(2)%Re(2)' : 0.0001, 'dt' : 1e-11,
            'viscous': 'T'})

        cases.append(define_case_d(stack, "",             {'weno_Re_flux': 'F'}))
        cases.append(define_case_d(stack, "weno_Re_flux", {'weno_Re_flux': 'T'}))
        for weno_Re_flux in ['T']:
            stack.push("weno_Re_flux" if weno_Re_flux == 'T' else '', {'weno_Re_flux' : 'T'})
            cases.append(define_case_d(stack, "weno_avg", {'weno_avg': 'T'}))
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

        cases.append(define_case_d(stack, "model_eqns=2", {'model_eqns': 2}))

        stack.push('cfl_adap_dt=T', {'cfl_adap_dt': 'T', 'cfl_target': 0.08, 't_save': 0.1, 'n_start': 0, 't_stop': 0.1})
        cases.append(define_case_d(stack, '', {}))

        stack.pop()

        stack.push("Viscous", {
            'fluid_pp(1)%Re(1)' : 0.0001, 'fluid_pp(1)%Re(2)' : 0.0001,
            'fluid_pp(2)%Re(1)' : 0.0001, 'fluid_pp(2)%Re(2)' : 0.0001, 'dt' : 1e-10,
            'viscous': 'T'
            })

        cases.append(define_case_d(stack, "",             {'weno_Re_flux': 'F'}))
        cases.append(define_case_d(stack, "weno_Re_flux", {'weno_Re_flux': 'T'}))
        for weno_Re_flux in ['T']:
            stack.push("weno_Re_flux" if weno_Re_flux == 'T' else '', {'weno_Re_flux' : 'T'})
            cases.append(define_case_d(stack, "weno_avg", {'weno_avg': 'T'}))
            stack.pop()

        stack.pop()
        stack.pop()

    def alter_ppn(dimInfo):
        if len(dimInfo[0]) == 3:
            cases.append(define_case_d(stack, '2 MPI Ranks', {'m': 29, 'n': 29, 'p': 49}, ppn=2))
            cases.append(define_case_d(stack, '2 MPI Ranks -> RDMA MPI', {'m': 29, 'n': 29, 'p': 49, 'rdma_mpi': 'T'}, ppn=2))
        else:
            cases.append(define_case_d(stack, '2 MPI Ranks', {}, ppn=2))
            cases.append(define_case_d(stack, '2 MPI Ranks -> RDMA MPI', {'rdma_mpi': 'T'}, ppn=2))

    def alter_ib(dimInfo, six_eqn_model=False):
        for slip in [True, False]:
            stack.push(f'IBM', {
                'ib': 'T', 'num_ibs': 1,
                'patch_ib(1)%x_centroid': 0.5, 'patch_ib(1)%y_centroid': 0.5,
                'patch_ib(1)%radius': 0.1, 'patch_icpp(1)%vel(1)': 0.001,
                'patch_icpp(2)%vel(1)': 0.001, 'patch_icpp(3)%vel(1)': 0.001,
                'patch_ib(1)%slip': 'T' if slip else 'F',
            })

            suffix = " -> slip" if slip else ""

            if len(dimInfo[0]) == 3:
                cases.append(define_case_d(stack, f'Sphere{suffix}', {
                    'patch_ib(1)%z_centroid': 0.5,
                    'patch_ib(1)%geometry': 8,
                }))

                cases.append(define_case_d(stack, f'Cuboid{suffix}', {
                    'patch_ib(1)%z_centroid': 0.5,
                    'patch_ib(1)%length_x': 0.1,
                    'patch_ib(1)%length_y': 0.1,
                    'patch_ib(1)%length_z': 0.1,
                    'patch_ib(1)%geometry': 9,
                }))

                cases.append(define_case_d(stack, f'Cylinder{suffix}', {
                    'patch_ib(1)%z_centroid': 0.5,
                    'patch_ib(1)%length_x': 0.1,
                    'patch_ib(1)%geometry': 10,
                }))

            elif len(dimInfo[0]) == 2:
                cases.append(define_case_d(stack, f'Rectangle{suffix}', {
                    'patch_ib(1)%length_x': 0.05,
                    'patch_ib(1)%length_y': 0.05,
                    'patch_ib(1)%geometry': 3,
                }))
                cases.append(define_case_d(stack, f'Circle{suffix}', {'patch_ib(1)%geometry': 2 }))
                if six_eqn_model:
                    cases.append(define_case_d(stack, f'model_eqns=3{suffix}', {
                        'patch_ib(1)%geometry': 2,
                        'model_eqns': 3,
                    }))

            stack.pop()

    def ibm_stl():
        common_mods = {
        't_step_stop': Nt, 't_step_save': Nt
        }
        for ndim in range(2, 4):
            cases.append(define_case_f(
                f'{ndim}D -> IBM -> STL',
                f'examples/{ndim}D_ibm_stl_test/case.py',
                ['--ndim', str(ndim)],
                mods=common_mods
            ))
    ibm_stl()

    def alter_acoustic_src(dimInfo):
        stack.push("Acoustic Source", {"acoustic_source": 'T', 'acoustic(1)%support': 1, 'dt': 1e-3, 't_step_stop': 50, 't_step_save': 50})

        transducer_params = {'acoustic(1)%loc(1)': 0.2, 'acoustic(1)%foc_length': 0.4, 'acoustic(1)%aperture': 0.6}

        if len(dimInfo[0]) == 1:
            for pulse_type in ['Sine', 'Square']:
                stack.push(pulse_type, {'acoustic(1)%pulse': 1 if pulse_type == 'Sine' else 3})
                cases.append(define_case_d(stack, 'Frequency', {'acoustic(1)%frequency': 50}))
                cases.append(define_case_d(stack, 'Wavelength', {'acoustic(1)%wavelength': 0.02}))
                cases.append(define_case_d(stack, 'Delay', {'acoustic(1)%delay': 0.02, 'acoustic(1)%wavelength': 0.02}))
                cases.append(define_case_d(stack, 'Number of Pulses', {'acoustic(1)%npulse': 2, 'acoustic(1)%wavelength': 0.01}))
                stack.pop()

            stack.push('Gaussian', {'acoustic(1)%pulse': 2, 'acoustic(1)%delay': 0.02})
            cases.append(define_case_d(stack, 'Sigma Time', {'acoustic(1)%gauss_sigma_time': 0.01}))
            cases.append(define_case_d(stack, 'Sigma Dist', {'acoustic(1)%gauss_sigma_dist': 0.01}))
            cases.append(define_case_d(stack, 'Dipole', {'acoustic(1)%gauss_sigma_dist': 0.01, 'acoustic(1)%dipole': 'T'}))
            stack.pop()

        elif len(dimInfo[0]) == 2:
            stack.push('', {'acoustic(1)%loc(2)': 0.5, 'acoustic(1)%wavelength': 0.02})

            stack.push('Planar', {})
            stack.push('support=2', {'acoustic(1)%support': 2})
            cases.append(define_case_d(stack, '', {}))
            cases.append(define_case_d(stack, 'Dipole', {'acoustic(1)%dipole': 'T'}))
            stack.pop()
            stack.pop()

            stack.push('Transducer', transducer_params)
            for support in [5, 6]:
                stack.push(f'support={support}', {'acoustic(1)%support': support, 'cyl_coord': 'T' if support == 6 else 'F', 'bc_y%beg': -2 if support == 6 else -3})
                cases.append(define_case_d(stack, 'Sine', {}))
                cases.append(define_case_d(stack, 'Gaussian', {'acoustic(1)%pulse': 2, 'acoustic(1)%delay': 0.02, 'acoustic(1)%gauss_sigma_dist': 0.01}))
                cases.append(define_case_d(stack, 'Delay', {'acoustic(1)%delay': 0.02}))
                stack.pop()
            stack.pop()

            stack.push('Transducer Array', {**transducer_params, 'acoustic(1)%num_elements': 4, 'acoustic(1)%element_spacing_angle': 0.05, 'acoustic(1)%element_on': 0})
            stack.push('support=9', {'acoustic(1)%support': 9})
            cases.append(define_case_d(stack, 'All Elements', {}))
            cases.append(define_case_d(stack, 'One element', {'acoustic(1)%element_on': 1}))
            stack.pop()
            cases.append(define_case_d(stack, 'support=10', {'acoustic(1)%support': 10, 'cyl_coord': 'T', 'bc_y%beg': -2}))
            stack.pop()

            stack.pop()

        elif len(dimInfo[0]) == 3:
            stack.push('', {'acoustic(1)%loc(2)': 0.5, 'acoustic(1)%loc(3)': 0.5, 'acoustic(1)%wavelength': 0.02})

            stack.push('Planar', {})
            stack.push('support=3', {'acoustic(1)%support': 3, 'acoustic(1)%height': 0.25})
            cases.append(define_case_d(stack, '', {}))
            cases.append(define_case_d(stack, 'Dipole', {'acoustic(1)%dipole': 'T'}))
            stack.pop()
            stack.pop()

            stack.push('Transducer', transducer_params)
            cases.append(define_case_d(stack, 'support=7', {'acoustic(1)%support': 7}))
            stack.pop()

            stack.push('Transducer Array', {**transducer_params, 'acoustic(1)%num_elements': 6, 'acoustic(1)%element_polygon_ratio': 0.7})
            stack.push('support=11', {'acoustic(1)%support': 11})
            cases.append(define_case_d(stack, 'All Elements', {}))
            cases.append(define_case_d(stack, 'One element', {'acoustic(1)%element_on': 1}))
            stack.pop()
            stack.pop()

            stack.pop()

        stack.pop()

    def alter_bubbles(dimInfo):
        if len(dimInfo[0]) > 0:
            stack.push("Bubbles", {"bubbles_euler": 'T'})

            stack.push('', {
                'nb' : 3, 'fluid_pp(1)%gamma' : 0.16, 'fluid_pp(1)%pi_inf': 3515.0,
                'fluid_pp(2)%gamma': 2.5, 'fluid_pp(2)%pi_inf': 0.0, 'fluid_pp(1)%mul0' : 0.001002,
                'fluid_pp(1)%ss' : 0.07275,'fluid_pp(1)%pv' : 2338.8,'fluid_pp(1)%gamma_v' : 1.33,
                'fluid_pp(1)%M_v' : 18.02,'fluid_pp(1)%mu_v' : 8.816e-06,'fluid_pp(1)%k_v' : 0.019426,
                'fluid_pp(2)%gamma_v' : 1.4,'fluid_pp(2)%M_v' : 28.97,'fluid_pp(2)%mu_v' : 1.8e-05,
                'fluid_pp(2)%k_v' : 0.02556, 'patch_icpp(1)%alpha_rho(1)': 0.96, 'patch_icpp(1)%alpha(1)':
                4e-02, 'patch_icpp(2)%alpha_rho(1)': 0.96, 'patch_icpp(2)%alpha(1)': 4e-02,  'patch_icpp(3)%alpha_rho(1)': 0.96,
                'patch_icpp(3)%alpha(1)': 4e-02, 'patch_icpp(1)%pres': 1.0, 'patch_icpp(2)%pres': 1.0,
                'patch_icpp(3)%pres': 1.0, 'acoustic(1)%support': 1, 'acoustic(1)%wavelength': 0.25
            })

            stack.push('', {"acoustic_source": 'T'})

            if len(dimInfo[0]) >= 2:
                stack.push("", {'acoustic(1)%loc(2)': 0.5, 'acoustic(1)%support': 2})

            if len(dimInfo[0]) >= 3:
                stack.push("", {'acoustic(1)%support': 3, 'acoustic(1)%height': 1e10})

            for polytropic in ['T', 'F']:
                stack.push("Polytropic" if polytropic == 'T' else '', {'polytropic' : polytropic})

                for bubble_model in [3, 2]:
                    stack.push(f"bubble_model={bubble_model}", {'bubble_model' : bubble_model})

                    if not (polytropic == 'F' and bubble_model == 3):
                        cases.append(define_case_d(stack, '', {}))

                    stack.pop()

                stack.pop()

            stack.push('', {'polytropic': 'T', 'bubble_model': 2})
            cases.append(define_case_d(stack, 'nb=1', {'nb': 1}))

            stack.push("adv_n=T", {'adv_n': 'T'})
            cases.append(define_case_d(stack, '', {}))
            cases.append(define_case_d(stack, 'adap_dt=T', {'adap_dt': 'T'}))
            stack.pop()

            stack.push('', {'fluid_pp(1)%pi_inf': 351.5})
            cases.append(define_case_d(stack, 'artificial_Ma', {'pi_fac': 0.1}))

            stack.pop()

            cases.append(define_case_d(stack, 'low_Mach=1', {'low_Mach': 1}))
            cases.append(define_case_d(stack, 'low_Mach=2', {'low_Mach': 2}))

            stack.push("QBMM", {'qbmm': 'T'})
            cases.append(define_case_d(stack, '', {}))


            stack.push("Non-polytropic", {'polytropic': 'F'})
            cases.append(define_case_d(stack, '', {}))

            stack.pop()

            stack.push('bubble_model=3', {'bubble_model': 3, 'polytropic': 'T'})
            cases.append(define_case_d(stack, '', {}))

            stack.push('Non-polytropic', { 'polytropic': 'F'})
            cases.append(define_case_d(stack, '', {}))

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
                'fd_order':                   4,
                'fluid_pp(1)%gamma':          0.3,    'fluid_pp(1)%pi_inf':         7.8E+05,
                'patch_icpp(1)%pres':         1.E+06, 'patch_icpp(1)%alpha_rho(1)': 1000.E+00,
                'patch_icpp(2)%pres':         1.E+05, 'patch_icpp(2)%alpha_rho(1)': 1000.E+00,
                'patch_icpp(3)%pres':         5.E+05, 'patch_icpp(3)%alpha_rho(1)': 1000.E+00,
                'patch_icpp(1)%tau_e(1)':     0.E-00, 'patch_icpp(2)%tau_e(1)':     0.E-00,
                'patch_icpp(3)%tau_e(1)':     0.E-00, 'fluid_pp(1)%G':              1.E+05,
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

            cases.append(define_case_d(stack, '', {}))

            reflective_params = {'bc_x%beg': -2, 'bc_x%end': -2, 'bc_y%beg': -2, 'bc_y%end': -2}
            if len(dimInfo[0]) == 3:
                reflective_params.update({'bc_z%beg': -2, 'bc_z%end': -2})

            if num_fluids == 1:
                cases.append(define_case_d(stack, 'cont_damage', {'cont_damage': 'T', 'tau_star': 0.0, 'cont_damage_s': 2.0, 'alpha_bar': 1e-4}))
                if len(dimInfo[0]) >= 2:
                    cases.append(define_case_d(stack, 'bc=-2', reflective_params))
                if len(dimInfo[0]) == 2:
                    cases.append(define_case_d(stack, 'Axisymmetric', {**reflective_params, 'cyl_coord': 'T'}))

            stack.pop()

            if num_fluids == 2:
                stack.pop()

            if len(dimInfo[0]) == 2:
                stack.pop()

            if len(dimInfo[0]) == 3:
                for _ in range(2):
                    stack.pop()

    def alter_body_forces(dimInfo):
        ndims = len(dimInfo[0])

        stack.push("Bodyforces",{
                'bf_x' : 'T', 'k_x' : 1, 'w_x' : 1, 'p_x': 1, 'g_x' : 10
            })

        if ndims >= 2:
            stack.push("",{
                    'bf_y' : 'T', 'k_y' : 1, 'w_y' : 1, 'p_y': 1, 'g_y' : 10
                })

        if ndims == 3:
            stack.push("",{
                    'bf_z' : 'T', 'k_z' : 1, 'w_z' : 1, 'p_z': 1, 'g_z' : 10
                })

        cases.append(define_case_d(stack, '', {}))

        stack.push('cfl_adap_dt=T', {'cfl_adap_dt': 'T', 'cfl_target': 0.08, 't_save': 0.025, 'n_start': 0, 't_stop': 0.025})
        cases.append(define_case_d(stack, '', {}))

        stack.pop()

        stack.pop()

        if ndims >= 2:
            stack.pop()

        if ndims == 3:
            stack.pop()

    def alter_mixlayer_perturb(dimInfo):
        if len(dimInfo[0]) == 3:
            cases.append(define_case_d(stack,'mixlayer_perturb',{
                'm': 24, 'n': 64, 'p': 24, 'dt': 1e-2,
                'num_patches': 1, 'num_fluids': 1,
                'x_domain%beg': 0.0, 'x_domain%end': 20.0, 'bc_x%beg': -1, 'bc_x%end': -1,
                'y_domain%beg': -10.0, 'y_domain%end': 10.0, 'bc_y%beg': -6, 'bc_y%end': -6,
                'z_domain%beg': 0.0, 'z_domain%end': 20.0, 'bc_z%beg': -1, 'bc_z%end': -1,
                'mixlayer_vel_profile': 'T', 'mixlayer_perturb': 'T',
                'weno_Re_flux': 'F', 'weno_avg': 'T', 'wenoz': 'T',
                'fluid_pp(1)%gamma': 2.5, 'fluid_pp(1)%pi_inf': 0.0,
                'fluid_pp(1)%Re(1)': 1.6881644098979287, 'viscous': 'T',
                'patch_icpp(1)%geometry': 9,
                'patch_icpp(1)%x_centroid': 10.0, 'patch_icpp(1)%length_x': 20.0,
                'patch_icpp(1)%y_centroid': 0.0, 'patch_icpp(1)%length_y': 20.0,
                'patch_icpp(1)%z_centroid': 10.0, 'patch_icpp(1)%length_z': 20.0,
                'patch_icpp(1)%vel(1)': 1.0, 'patch_icpp(1)%vel(2)': 0.0, 'patch_icpp(1)%vel(3)': 0.0,
                'patch_icpp(1)%pres': 17.8571428571, 'patch_icpp(1)%alpha_rho(1)': 1.0, 'patch_icpp(1)%alpha(1)': 1.0,
                'patch_icpp(1)%r0': -1e6, 'patch_icpp(1)%v0': -1e6,
                'patch_icpp(2)%geometry': -100, 
                'patch_icpp(2)%x_centroid': -1e6, 'patch_icpp(2)%length_x': -1e6,
                'patch_icpp(2)%y_centroid': -1e6, 'patch_icpp(2)%length_y': -1e6, 
                'patch_icpp(2)%z_centroid': -1e6, 'patch_icpp(2)%length_z': -1e6, 
                'patch_icpp(2)%vel(1)': -1e6, 'patch_icpp(2)%vel(2)': -1e6, 'patch_icpp(2)%vel(3)': -1e6, 
                'patch_icpp(2)%r0': -1e6, 'patch_icpp(2)%v0': -1e6,
                'patch_icpp(3)%geometry': -100, 
                'patch_icpp(3)%x_centroid': -1e6, 'patch_icpp(3)%length_x': -1e6,
                'patch_icpp(3)%y_centroid': -1e6, 'patch_icpp(3)%length_y': -1e6, 
                'patch_icpp(3)%z_centroid': -1e6, 'patch_icpp(3)%length_z': -1e6, 
                'patch_icpp(3)%vel(1)': -1e6, 'patch_icpp(3)%vel(2)': -1e6, 'patch_icpp(3)%vel(3)': -1e6, 
                'patch_icpp(3)%r0': -1e6, 'patch_icpp(3)%v0': -1e6
            }))

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

                cases.append(define_case_d(stack, '', {}))

                stack.pop()
                stack.pop()

                if num_fluids == 3:
                    stack.pop()

    def alter_viscosity(dimInfo):
        # Viscosity & bubbles checks
        if len(dimInfo[0]) > 0:
            stack.push("Viscosity -> Bubbles",
                       {"fluid_pp(1)%Re(1)": 50, "bubbles_euler": 'T', "viscous": 'T'})

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
                        cases.append(define_case_d(stack, '', {}))

                    stack.pop()

                stack.pop()

            stack.push('', {'polytropic': 'T', 'bubble_model': 2})
            cases.append(define_case_d(stack, 'nb=1', {'nb': 1}))

            stack.push("QBMM", {'qbmm': 'T'})
            cases.append(define_case_d(stack, '', {}))

            stack.push('bubble_model=3', {'bubble_model': 3})
            cases.append(define_case_d(stack, '', {}))

            stack.push('cfl_adap_dt=T', {'cfl_adap_dt': 'T', 'cfl_target': 0.8, 't_save': 0.01, 'n_start': 0, 't_stop': 0.01, 'm': 24})
            cases.append(define_case_d(stack, '', {}))

            stack.pop()

            stack.push('cfl_const_dt=T', {'cfl_const_dt': 'T', 'cfl_target': 0.8, 't_save': 0.01, 'n_start': 0, 't_stop': 0.01, 'm': 24})
            cases.append(define_case_d(stack, '', {}))

            for _ in range(6):
                stack.pop()

    def alter_lag_bubbles():
        # Lagrangian bubbles
        for adap_dt in ['F', 'T']:
            for couplingMethod in [1, 2]:
                stack.push("Lagrange Bubbles", {"bubbles_lagrange": 'T',
                    'dt': 1e-06, 'lag_params%pressure_corrector': 'T', 'bubble_model': 2,
                    'num_fluids': 2, 'lag_params%heatTransfer_model': 'T', 'lag_params%massTransfer_model': 'T',
                    'fluid_pp(1)%gamma' : 0.16, 'fluid_pp(1)%pi_inf': 3515.0, 'fluid_pp(2)%gamma': 2.5,
                    'fluid_pp(2)%pi_inf': 0.0, 'fluid_pp(1)%mul0' : 0.001002, 'fluid_pp(1)%ss' : 0.07275,
                    'fluid_pp(1)%pv' : 2338.8,'fluid_pp(1)%gamma_v' : 1.33, 'fluid_pp(1)%M_v' : 18.02,
                    'fluid_pp(1)%mu_v' : 8.816e-06,'fluid_pp(1)%k_v' : 0.019426, 'fluid_pp(1)%cp_v' : 2.1e3,
                    'fluid_pp(2)%gamma_v' : 1.4,'fluid_pp(2)%M_v' : 28.97, 'fluid_pp(2)%mu_v' : 1.8e-05,
                    'fluid_pp(2)%k_v' : 0.02556, 'fluid_pp(2)%cp_v' : 1.e3, 'patch_icpp(1)%alpha_rho(1)': 0.96,
                    'patch_icpp(1)%alpha(1)': 4e-02, 'patch_icpp(1)%alpha_rho(2)': 0., 'patch_icpp(1)%alpha(2)': 0.,
                    'patch_icpp(2)%alpha_rho(1)': 0.96, 'patch_icpp(2)%alpha(1)': 4e-02, 'patch_icpp(2)%alpha_rho(2)': 0.,
                    'patch_icpp(2)%alpha(2)': 0.,  'patch_icpp(3)%alpha_rho(1)': 0.96, 'patch_icpp(3)%alpha(1)': 4e-02,
                    'patch_icpp(3)%alpha_rho(2)': 0., 'patch_icpp(3)%alpha(2)': 0.,'patch_icpp(1)%pres': 1.0,
                    'patch_icpp(2)%pres': 1.0, 'patch_icpp(3)%pres': 1.0, 'acoustic_source': 'T', 'acoustic(1)%loc(2)': 0.5,
                    'acoustic(1)%wavelength': 0.25, 'acoustic(1)%support': 3, 'acoustic(1)%height': 1e10,
                    'acoustic(1)%mag': 2e+04, 't_step_start': 0, 't_step_stop': 50, 't_step_save': 50
                })
                if couplingMethod==1:
                    stack.push('One-way Coupling',{'lag_params%solver_approach': 1})
                else:
                    stack.push('Two-way Coupling',{'lag_params%solver_approach': 2})

                if adap_dt=='F':
                    stack.push('',{})
                else:
                    stack.push('adap_dt=T',{'adap_dt': 'T'})

                cases.append(define_case_d(stack, '', {}))

                stack.pop()

                stack.pop()

                stack.pop()

    def alter_elliptic_smoothing():
        # Elliptic Smoothing

        stack.push("Smoothing",{
                'elliptic_smoothing': 'T', 'elliptic_smoothing_iters': 10
            })

        cases.append(define_case_d(stack, '', {}))

        stack.pop()

    def alter_bc_patches(dimInfo):
       # BC_Patches

        stack.push('BC Patches',{
            'num_bc_patches': 1
        })

        if len(dimInfo[0]) > 2:
            for direc in [1,2,3]:
                stack.push('Circle' ,{
                        'patch_bc(1)%geometry': 2, 'patch_bc(1)%dir': direc,
                        'patch_bc(1)%type': -17, 'patch_bc(1)%loc': -1
                })

                if direc==1:
                    stack.push('X', {'patch_bc(1)%centroid(2)': 0,'patch_bc(1)%centroid(3)': 0, "patch_bc(1)%radius": 0.000125,})
                elif direc==2:
                    stack.push('Y', {'patch_bc(1)%centroid(1)': 0,'patch_bc(1)%centroid(3)': 0, "patch_bc(1)%radius": 0.000125,})
                else:
                    stack.push('Z', {'patch_bc(1)%centroid(1)': 0,'patch_bc(1)%centroid(2)': 0, "patch_bc(1)%radius": 0.000125,})


                cases.append(define_case_d(stack, '', {}))

                stack.pop()

                stack.pop()

        elif len(dimInfo[0]) > 1:
            for direc in [1,2]:
                stack.push('Line Segment' ,{
                        'patch_bc(1)%geometry': 1, 'patch_bc(1)%dir': direc,
                        'patch_bc(1)%type': -17, 'patch_bc(1)%loc': -1
                })

                if direc==1:
                    stack.push('X' ,{'patch_bc(1)%centroid(2)': 0.0, 'patch_bc(1)%length(2)': 0.0025})
                else:
                    stack.push('Y' ,{'patch_bc(1)%centroid(1)': 0.0, 'patch_bc(1)%length(1)': 0.0025})


                cases.append(define_case_d(stack, '', {}))

                stack.pop()

                stack.pop()

        stack.pop()

    def mhd_cases():
        params = {
            '1D': {"m": 200, "dt": 0.001, "t_step_stop": 200, "t_step_save": 200},
            '2D': {"m": 50, "n": 50, "dt": 0.002, "t_step_stop": 500, "t_step_save": 500},
            '3D': {"m": 25, "n": 25, "p": 25, "dt": 0.005, "t_step_stop": 200, "t_step_save": 200},
        }

        case_specs = [
            ("1D -> MHD -> HLL",    "examples/1D_brio_wu/case.py",              params['1D']),
            ("1D -> MHD -> HLLD",   "examples/1D_brio_wu_hlld/case.py",         params['1D']),
            ("1D -> RMHD",          "examples/1D_brio_wu_rmhd/case.py",         params['1D']),
            ("2D -> MHD -> HLL",    "examples/2D_orszag_tang/case.py",          params['2D']),
            ("2D -> MHD -> HLLD",   "examples/2D_orszag_tang/case.py",          {**params['2D'], 'riemann_solver': 4}),
            ("2D -> MHD -> Powell", "examples/2D_orszag_tang_powell/case.py",   params['2D']),
            ("2D -> RMHD",          "examples/2D_shock_cloud_rmhd/case.py",     params['2D']),
            ("3D -> MHD",           "examples/3D_brio_wu/case.py",              params['3D']),
            ("3D -> RMHD",          "examples/3D_brio_wu/case.py",              {**params['3D'], 'relativity': 'T'}),
        ]

        for name, path, param in case_specs:
            cases.append(define_case_f(name, path, mods=param))

    def foreach_dimension():
        for dimInfo, dimParams in get_dimensions():
            stack.push(f"{len(dimInfo[0])}D", dimParams)
            alter_bcs(dimInfo)
            alter_grcbc(dimInfo)
            alter_weno(dimInfo)
            alter_muscl()
            alter_num_fluids(dimInfo)
            if len(dimInfo[0]) == 2:
                alter_2d()
            if len(dimInfo[0]) == 3:
                alter_3d()
                alter_lag_bubbles()
            alter_ppn(dimInfo)
            stack.push('', {'dt': [1e-07, 1e-06, 1e-06][len(dimInfo[0])-1]})
            alter_acoustic_src(dimInfo)
            alter_bubbles(dimInfo)
            alter_hypoelasticity(dimInfo)
            alter_phasechange(dimInfo)
            alter_viscosity(dimInfo)
            alter_elliptic_smoothing()
            alter_body_forces(dimInfo)
            alter_mixlayer_perturb(dimInfo)
            alter_bc_patches(dimInfo)
            stack.pop()
            stack.pop()

    def foreach_example():
        for path in os.listdir(common.MFC_EXAMPLE_DIRPATH):
            if path == "scaling":
                continue

            # # List of all example cases that will be skipped during testing
            casesToSkip = ["2D_ibm_cfl_dt", "1D_sodHypo", "2D_viscous",
                           "2D_laplace_pressure_jump", "2D_bubbly_steady_shock",
                           "2D_advection", "2D_hardcoded_ic",
                           "2D_ibm_multiphase", "2D_acoustic_broadband",
                           "1D_inert_shocktube", "1D_reactive_shocktube",
                           "2D_ibm_steady_shock", "3D_performance_test",
                           "3D_ibm_stl_ellipsoid", "3D_sphbubcollapse",
                           "2D_ibm_stl_wedge", "3D_ibm_stl_pyramid",
                           "3D_ibm_bowshock", "3D_turb_mixing",
                           "2D_mixing_artificial_Ma",
                           "2D_lagrange_bubblescreen",
                           "3D_lagrange_bubblescreen", "2D_triple_point",
                           "1D_shuosher_analytical",
                           "1D_titarevtorro_analytical", 
                           "2D_acoustic_pulse_analytical",
                           "2D_isentropicvortex_analytical",
                           "2D_zero_circ_vortex_analytical",
                           "3D_TaylorGreenVortex_analytical",
                           "2D_backward_facing_step",
                           "2D_forward_facing_step"]
            if path in casesToSkip:
                continue
            name = f"{path.split('_')[0]} -> Example -> {'_'.join(path.split('_')[1:])}"
            path = os.path.join(common.MFC_EXAMPLE_DIRPATH, path, "case.py")
            if not os.path.isfile(path):
                continue
            def modify_example_case(case: dict):
                case['parallel_io'] = 'F'
                if 't_step_stop' in case and case['t_step_stop'] >= 50:
                    case['t_step_start'] = 0
                    case['t_step_stop'] = 50
                    case['t_step_save'] = 50

                caseSize = case['m'] * max(case['n'], 1) * max(case['p'], 1)
                if caseSize > 25 * 25:
                    if case['n'] == 0 and case['p'] == 0:
                        case['m'] = 25 * 25
                    elif case['p'] == 0:
                        case['m'] = 25
                        case['n'] = 25
                    elif caseSize > 25 * 25 * 25:
                        case['m'] = 25
                        case['n'] = 25
                        case['p'] = 25

            cases.append(define_case_f(name, path, [], {}, functor=modify_example_case))

    def chemistry_cases():
        common_mods = {
            't_step_stop': Nt, 't_step_save': Nt
        }
        for ndim in range(1, 4):
            cases.append(define_case_f(
                f'{ndim}D -> Chemistry -> Perfect Reactor',
                'examples/nD_perfect_reactor/case.py',
                ['--ndim', str(ndim)],
                mods=common_mods
            ))

        for riemann_solver, gamma_method in itertools.product([1, 2], [1, 2]):
            cases.append(define_case_f(
                f'1D -> Chemistry -> Inert Shocktube -> Riemann Solver {riemann_solver} -> Gamma Method {gamma_method}',
                'examples/1D_inert_shocktube/case.py',
                mods={
                    **common_mods,
                    'riemann_solver': riemann_solver,
                    'chem_params%gamma_method': gamma_method
                },
                override_tol=1
            ))

    foreach_dimension()

    mhd_cases()

    foreach_example()

    chemistry_cases()

    # Sanity Check 1
    if stack.size() != 0:
        raise common.MFCException("list_cases: stack isn't fully pop'ed")

    # Sanity Check 2
    uuids  = [ case.get_uuid() for case in cases ]
    l1, l2 = len(uuids), len(set(uuids))
    if l1 != l2:
        raise common.MFCException(f"list_cases: uuids aren't unique ({l1} cases but {l2} unique uuids)")

    return cases
