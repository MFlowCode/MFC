import dataclasses
import glob
import json
import os
import typing

from .. import case_validator, common
from ..case import Case

# Note: pyrometheus and cantera are imported lazily in the methods that need them
# to avoid slow startup times for commands that don't use chemistry features
# Note: build is imported lazily to avoid circular import with build.py
from ..printer import cons
from ..state import ARG, ARGS, gpuConfigOptions


@dataclasses.dataclass(init=False)
class MFCInputFile(Case):
    filename: str
    dirpath: str

    def __init__(self, filename: str, dirpath: str, params: dict) -> None:
        super().__init__(params)
        self.filename = filename
        self.dirpath = dirpath

    def generate_inp(self, target) -> None:
        from .. import build

        target = build.get_target(target)

        # Save .inp input file
        common.file_write(f"{self.dirpath}/{target.name}.inp", self.get_inp(target))

    def __save_fpp(self, target, contents: str) -> None:
        inc_dir = os.path.join(target.get_staging_dirpath(self), "include", target.name)
        common.create_directory(inc_dir)

        fpp_path = os.path.join(inc_dir, "case.fpp")

        cons.print("Writing a (new) custom case.fpp file.")
        common.file_write(fpp_path, contents, True)

    def get_cantera_solution(self):
        # Lazy import to avoid slow startup for commands that don't need chemistry
        import cantera as ct

        if self.params.get("chemistry", "F") == "T":
            cantera_file = self.params["cantera_file"]
            candidates = [
                cantera_file,
                os.path.join(self.dirpath, cantera_file),
                os.path.join(common.MFC_MECHANISMS_DIR, cantera_file),
            ]
        else:
            # Chemistry is off — return a dummy solution so MFC still compiles.
            cantera_file = "h2o2.yaml"
            candidates = [cantera_file]

        for candidate in candidates:
            try:
                return ct.Solution(candidate)
            except Exception as e:
                cons.print(f"[dim]  Cantera: skipping '{candidate}': {e}[/dim]")
                continue

        raise common.MFCException(f"Cantera file '{cantera_file}' not found. Searched: {', '.join(candidates)}.")

    def get_cantera_surface(self):
        # Lazy import to avoid slow startup for commands that don't need chemistry
        import cantera as ct
        import yaml

        surface_file = self.params.get("surface_cantera_file")
        surface_phase = self.params.get("surface_phase")

        if surface_file is None and surface_phase is None:
            return None

        if surface_file is None or surface_phase is None:
            raise common.MFCException("surface_cantera_file and surface_phase must be specified together.")

        candidates = [
            surface_file,
            os.path.join(self.dirpath, surface_file),
            os.path.join(common.MFC_MECHANISMS_DIR, surface_file),
        ]

        gas = self.get_cantera_solution()

        for candidate in candidates:
            if not os.path.isfile(candidate):
                continue

            try:
                with open(candidate, "r", encoding="utf-8") as stream:
                    mechanism = yaml.safe_load(stream)

                phases = mechanism.get("phases", [])

                interface_data = None
                for phase in phases:
                    if phase.get("name") == surface_phase:
                        interface_data = phase
                        break

                if interface_data is None:
                    raise common.MFCException(f"Surface phase '{surface_phase}' was not found in '{candidate}'.")

                adjacent_names = interface_data.get("adjacent-phases", [])

                adjacent = []

                for phase_name in adjacent_names:
                    if phase_name == gas.name:
                        adjacent.append(gas)
                    else:
                        adjacent.append(ct.Solution(candidate, phase_name))

                return ct.Interface(
                    candidate,
                    surface_phase,
                    adjacent=adjacent,
                )

            except common.MFCException:
                raise
            except Exception as e:
                cons.print(f"[dim]  Cantera: skipping surface mechanism " f"'{candidate}': {e}[/dim]")
                continue

        raise common.MFCException(f"Cantera surface file '{surface_file}' with phase " f"'{surface_phase}' could not be loaded. " f"Searched: {', '.join(candidates)}.")

    def generate_surface_thermochem(self, sol, surface, directive_str=None) -> str:
        """Generate the MFC heterogeneous surface-chemistry Fortran module."""

        if directive_str == "mp":
            gpu_routine_define = "#define GPU_ROUTINE(name) !$omp declare target device_type(any)"
        elif directive_str == "acc":
            gpu_routine_define = "#define GPU_ROUTINE(name) !$acc routine seq"
        else:
            gpu_routine_define = "#define GPU_ROUTINE(name) ! name"

        if surface is None:
            num_species = len(sol.species_names)

            lines = [
                "! This file is automatically generated by the MFC toolchain.",
                "! Do not edit manually.",
                "",
                gpu_routine_define,
                "",
                "module m_surface_thermochem",
                "",
                "    use m_precision_select, only: wp",
                "",
                "    implicit none",
                "",
                "    private",
                "    public :: get_surface_net_production_rates",
                "    public :: get_surface_reaction_heat_flux",
                "",
                "contains",
                "",
                "    subroutine get_surface_net_production_rates( &",
                "        density, temperature, mass_fractions, omega_s)",
                "",
                "        GPU_ROUTINE(get_surface_net_production_rates)",
                "",
                "        real(wp), intent(in) :: density",
                "        real(wp), intent(in) :: temperature",
                f"        real(wp), intent(in) :: mass_fractions({num_species})",
                f"        real(wp), intent(out) :: omega_s({num_species})",
                "",
                "        omega_s = 0._wp",
                "",
                "    end subroutine get_surface_net_production_rates",
                "",
                "    subroutine get_surface_reaction_heat_flux( &",
                "        density, temperature, mass_fractions, q_rxn)",
                "",
                "        GPU_ROUTINE(get_surface_reaction_heat_flux)",
                "",
                "        real(wp), intent(in) :: density",
                "        real(wp), intent(in) :: temperature",
                f"        real(wp), intent(in) :: mass_fractions({num_species})",
                "        real(wp), intent(out) :: q_rxn",
                "",
                "        q_rxn = 0._wp",
                "",
                "    end subroutine get_surface_reaction_heat_flux",
                "",
                "end module m_surface_thermochem",
                "",
            ]

            return "\n".join(lines)

        gas_species = sol.species_names
        gas_index = {name: i + 1 for i, name in enumerate(gas_species)}

        # Collect non-gas species participating in heterogeneous reactions.
        # Their thermodynamic data are needed for reaction enthalpies.
        nongas_species = {}

        for phase in surface.adjacent.values():
            for species in phase.species():
                if species.name not in gas_index:
                    nongas_species[species.name] = species

        # Verify that all species appearing in the surface reactions belong to
        # either the gas mechanism, the surface phase, or an adjacent phase.
        for reaction in surface.reactions():
            for name in set(reaction.reactants) | set(reaction.products):
                if name in gas_index:
                    continue
                if name in surface.species_names:
                    raise common.MFCException(f"Surface-site species '{name}' in heterogeneous reaction " "stoichiometry are not currently supported.")
                if name in nongas_species:
                    continue

                raise common.MFCException(f"Surface reaction species '{name}' is not present in any " "phase associated with the surface mechanism.")

        # Return a Fortran expression for h/(R*T) for a non-gas species.
        # Cantera's NasaPoly2 stores:
        #   [Tmid, a1_high, ..., a7_high, a1_low, ..., a7_low]
        # and
        #   h/(R*T) = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4
        #             + a5*T^4/5 + a6/T.
        def nongas_h_rt_expressions(species):
            thermo = species.thermo

            if thermo is None or thermo.__class__.__name__ != "NasaPoly2":
                raise common.MFCException(f"Surface thermochemistry for non-gas species " f"'{species.name}' requires Cantera NasaPoly2 thermo data.")

            coeffs = list(thermo.coeffs)
            if len(coeffs) != 15:
                raise common.MFCException(f"Unexpected NasaPoly2 coefficient count for " f"surface species '{species.name}'.")

            tmid = coeffs[0]
            high = coeffs[1:8]
            low = coeffs[8:15]

            def h_rt(a):
                return (
                    f"({a[0]:.16e}_wp"
                    f" + ({a[1]:.16e}_wp)*temperature/2._wp"
                    f" + ({a[2]:.16e}_wp)*temperature**2/3._wp"
                    f" + ({a[3]:.16e}_wp)*temperature**3/4._wp"
                    f" + ({a[4]:.16e}_wp)*temperature**4/5._wp"
                    f" + ({a[5]:.16e}_wp)/temperature)"
                )

            return tmid, h_rt(low), h_rt(high)

        # Precompute generated h/(R*T) expressions for participating
        # non-gas species from adjacent phases.
        nongas_h_rt = {}
        reaction_species = set()
        for reaction in surface.reactions():
            reaction_species.update(reaction.reactants)
            reaction_species.update(reaction.products)

        for name in reaction_species:
            if name in gas_index:
                continue
            if name in nongas_species:
                nongas_h_rt[name] = nongas_h_rt_expressions(nongas_species[name])

        lines = [
            "! This file is automatically generated by the MFC toolchain.",
            "! Do not edit manually.",
            "",
            gpu_routine_define,
            "",
            "module m_surface_thermochem",
            "",
            "    use m_precision_select, only: wp",
            "    use m_thermochem, only: gas_constant, get_species_enthalpies_rt",
            "",
            "    implicit none",
            "",
            "    private",
            "    public :: get_surface_net_production_rates",
            "    public :: get_surface_reaction_heat_flux",
            "",
            "contains",
            "",
            "    subroutine get_surface_net_production_rates( &",
            "        density, temperature, mass_fractions, omega_s)",
            "",
            "        GPU_ROUTINE(get_surface_net_production_rates)",
            "",
            "        real(wp), intent(in) :: density",
            "        real(wp), intent(in) :: temperature",
            f"        real(wp), intent(in) :: mass_fractions({len(gas_species)})",
            f"        real(wp), intent(out) :: omega_s({len(gas_species)})",
            "",
            f"        real(wp) :: concentrations({len(gas_species)})",
            "        real(wp) :: rate",
            "",
            "        omega_s = 0._wp",
            "        concentrations = 0._wp",
        ]

        for i, species in enumerate(sol.species()):
            lines.append(f"        concentrations({i + 1}) = " f"density*mass_fractions({i + 1})/" f"{species.molecular_weight:.16e}_wp")

        lines.append("")

        def append_reaction_rate(lines, reaction_number, reaction):
            rate = reaction.rate

            if not hasattr(rate, "pre_exponential_factor"):
                raise common.MFCException(f"Surface reaction {reaction_number} does not use a " "supported Arrhenius rate expression.")

            A = rate.pre_exponential_factor
            b = rate.temperature_exponent
            Ea = rate.activation_energy

            lines.extend(
                [
                    f"        ! Surface reaction {reaction_number}",
                    f"        rate = {A:.16e}_wp",
                ]
            )

            if b != 0.0:
                lines.append(f"        rate = rate*temperature**({b:.16e}_wp)")

            if Ea != 0.0:
                lines.append(f"        rate = rate*exp(-{Ea:.16e}_wp/" "(gas_constant*temperature))")

            # Cantera reaction orders override stoichiometric reactant orders.
            orders = dict(reaction.orders)
            for name in orders:
                if name in surface.species_names:
                    raise common.MFCException(f"Surface reaction {reaction_number} specifies an explicit " f"reaction order for surface-site species '{name}', which is " "not currently supported.")

            for name, nu in reaction.reactants.items():
                if name not in gas_index:
                    continue

                order = orders.get(name, nu)
                lines.append(f"        rate = rate*concentrations({gas_index[name]})**" f"({order:.16e}_wp)")

            # Explicit orders may include gas species not appearing as reactants.
            for name, order in orders.items():
                if name not in gas_index or name in reaction.reactants:
                    continue

                lines.append(f"        rate = rate*concentrations({gas_index[name]})**" f"({order:.16e}_wp)")

        for reaction_number, reaction in enumerate(surface.reactions(), start=1):
            append_reaction_rate(lines, reaction_number, reaction)

            for name in gas_species:
                nu = reaction.products.get(name, 0.0) - reaction.reactants.get(name, 0.0)
                if nu == 0.0:
                    continue

                lines.append(f"        omega_s({gas_index[name]}) = " f"omega_s({gas_index[name]}) " f"+ ({nu:.16e}_wp)*rate")

            lines.append("")

        lines.extend(
            [
                "    end subroutine get_surface_net_production_rates",
                "",
                "    subroutine get_surface_reaction_heat_flux( &",
                "        density, temperature, mass_fractions, q_rxn)",
                "",
                "        GPU_ROUTINE(get_surface_reaction_heat_flux)",
                "",
                "        real(wp), intent(in) :: density",
                "        real(wp), intent(in) :: temperature",
                f"        real(wp), intent(in) :: mass_fractions({len(gas_species)})",
                "        real(wp), intent(out) :: q_rxn",
                "",
                f"        real(wp) :: h0_rt({len(gas_species)})",
                f"        real(wp) :: concentrations({len(gas_species)})",
                "        real(wp) :: rate",
                "        real(wp) :: delta_h",
                "        real(wp) :: h_rt_nongas",
                "",
                "        call get_species_enthalpies_rt(temperature, h0_rt)",
                "",
                "        concentrations = 0._wp",
            ]
        )

        for i, species in enumerate(sol.species()):
            lines.append(f"        concentrations({i + 1}) = " f"density*mass_fractions({i + 1})/" f"{species.molecular_weight:.16e}_wp")

        lines.extend(
            [
                "",
                "        q_rxn = 0._wp",
                "",
            ]
        )

        for reaction_number, reaction in enumerate(surface.reactions(), start=1):
            append_reaction_rate(lines, reaction_number, reaction)
            lines.append("        delta_h = 0._wp")

            # Gas species contribution to reaction enthalpy:
            # h_k = (h_k/RT) * R * T, with R in J/(kmol K).
            for name in gas_species:
                nu = reaction.products.get(name, 0.0) - reaction.reactants.get(name, 0.0)
                if nu == 0.0:
                    continue

                lines.append(f"        delta_h = delta_h + ({nu:.16e}_wp)" f"*h0_rt({gas_index[name]})" "*gas_constant*temperature")

            # Non-gas species contribution. For Kestel this is graphite.
            for name in nongas_species:
                nu = reaction.products.get(name, 0.0) - reaction.reactants.get(name, 0.0)
                if nu == 0.0:
                    continue

                if name not in nongas_h_rt:
                    raise common.MFCException(f"No supported thermochemistry is available for " f"non-gas surface-reaction species '{name}'.")

                tmid, low_expr, high_expr = nongas_h_rt[name]

                lines.extend(
                    [
                        f"        if (temperature <= {tmid:.16e}_wp) then",
                        f"            h_rt_nongas = {low_expr}",
                        "        else",
                        f"            h_rt_nongas = {high_expr}",
                        "        end if",
                        f"        delta_h = delta_h + ({nu:.16e}_wp)" "*h_rt_nongas*gas_constant*temperature",
                    ]
                )

            # Positive q_rxn denotes heat released by an exothermic reaction.
            lines.append("        q_rxn = q_rxn - delta_h*rate")
            lines.append("")

        lines.extend(
            [
                "    end subroutine get_surface_reaction_heat_flux",
                "",
                "end module m_surface_thermochem",
                "",
            ]
        )

        return "\n".join(lines)

    def generate_fpp(self, target) -> None:
        # Lazy import to avoid slow startup for commands that don't need chemistry
        import pyrometheus as pyro

        if target.isDependency:
            return

        cons.print("Generating [magenta]case.fpp[/magenta].")
        cons.indent()

        # Case FPP file
        self.__save_fpp(target, self.get_fpp(target))

        # (Thermo)Chemistry source file
        modules_dir = os.path.join(target.get_staging_dirpath(self), "modules", target.name)
        common.create_directory(modules_dir)

        # Determine the real type based on the single precision flag
        real_type = "real(sp)" if (ARG("single") or ARG("mixed")) else "real(dp)"

        if ARG("gpu") == gpuConfigOptions.MP.value:
            directive_str = "mp"
        elif ARG("gpu") == gpuConfigOptions.ACC.value:
            directive_str = "acc"
        else:
            directive_str = None

        # Write the generated Fortran code to the m_thermochem.f90 file with the chosen precision
        sol = self.get_cantera_solution()
        surface = self.get_cantera_surface() if target.name == "simulation" else None
        if surface is not None:
            cons.print(f"Loaded Cantera surface phase '{surface.name}' " f"with {surface.n_reactions} reaction(s).")

        thermochem_code = pyro.FortranCodeGenerator().generate("m_thermochem", sol, pyro.CodeGenerationOptions(scalar_type=real_type, directive_offload=directive_str))

        common.file_write(os.path.join(modules_dir, "m_thermochem.f90"), thermochem_code, True)

        if target.name == "simulation":
            surface_thermochem_code = self.generate_surface_thermochem(
                sol,
                surface,
                directive_str,
            )

            surface_thermochem_path = os.path.join(
                modules_dir,
                "m_surface_thermochem.f90",
            )

            common.file_write(
                surface_thermochem_path,
                surface_thermochem_code,
                True,
            )

            if surface is not None:
                cons.print(f"Generated m_surface_thermochem.f90 with " f"{surface.n_reactions} surface reaction(s).")

        cons.unindent()

    def validate_constraints(self, target) -> None:
        """Validate case parameter constraints for a given target stage"""
        from .. import build

        target_obj = build.get_target(target)
        stage = target_obj.name

        try:
            warnings = case_validator.validate_case_constraints(self.params, stage)
        except case_validator.CaseConstraintError as e:
            raise common.MFCException(f"Case validation failed for {stage}:\n{e}") from e

        if warnings:
            cons.print()
            cons.print("[bold yellow]Physics warnings:[/bold yellow]")
            for warning in warnings:
                cons.print(f"  [yellow]- {warning}[/yellow]")
            cons.print()

    # Generate case.fpp & [target.name].inp
    def generate(self, target) -> None:
        # Validate constraints before generating input files
        self.validate_constraints(target)
        self.generate_inp(target)
        cons.print()
        self.generate_fpp(target)

    def clean(self, _targets) -> None:
        from .. import build

        targets = [build.get_target(target) for target in _targets]

        files = set()
        dirs = set()

        files = set(["equations.dat", "run_time.inf", "time_data.dat", "io_time_data.dat", "fort.1", "pre_time_data.dat"] + [f"{target.name}.inp" for target in targets])

        if build.PRE_PROCESS in targets:
            files = files | set(glob.glob(os.path.join(self.dirpath, "D", "*.000000.dat")))
            dirs = dirs | set(glob.glob(os.path.join(self.dirpath, "p_all", "p*", "0")))

        if build.SIMULATION in targets:
            restarts = set(glob.glob(os.path.join(self.dirpath, "restart_data", "*.dat")))
            restarts = restarts - set(glob.glob(os.path.join(self.dirpath, "restart_data", "lustre_0.dat")))
            restarts = restarts - set(glob.glob(os.path.join(self.dirpath, "restart_data", "lustre_*_cb.dat")))

            Ds = set(glob.glob(os.path.join(self.dirpath, "D", "*.dat")))
            Ds = Ds - set(glob.glob(os.path.join(self.dirpath, "D", "*.000000.dat")))

            files = files | restarts
            files = files | Ds

        if build.POST_PROCESS in targets:
            dirs.add("silo_hdf5")

        for relfile in files:
            filepath = relfile if os.path.isfile(relfile) else os.path.join(self.dirpath, relfile)
            common.delete_file(filepath)

        for reldir in dirs:
            dirpath = reldir if os.path.isdir(reldir) else os.path.join(self.dirpath, reldir)
            common.delete_directory(dirpath)


# Load the input file
def load(filepath: str = None, args: typing.List[str] = None, empty_data: dict = None, do_print: bool = True) -> MFCInputFile:
    if not filepath:
        if empty_data is None:
            raise common.MFCException("Please provide an input file.")

        input_file = MFCInputFile("empty.py", "empty.py", empty_data)
        input_file.validate_params()
        return input_file

    filename: str = filepath.strip()

    if do_print:
        cons.print(f"Acquiring [bold magenta]{filename}[/bold magenta]...")

    dirpath: str = os.path.abspath(os.path.dirname(filename))
    dictionary: dict = {}

    if not os.path.exists(filename):
        raise common.MFCException(f"Input file '{filename}' does not exist. Please check the path is valid.")

    if filename.endswith(".py"):
        (json_str, err) = common.get_py_program_output(filename, ["--mfc", json.dumps(ARGS())] + (args or []))

        if err != 0:
            raise common.MFCException(f"Input file {filename} terminated with a non-zero exit code. Please make sure running the file doesn't produce any errors.")
    elif filename.endswith(".json"):
        json_str = common.file_read(filename)
    elif filename.endswith((".yaml", ".yml")):
        import yaml

        with open(filename, "r") as f:
            dictionary = yaml.safe_load(f)
        json_str = json.dumps(dictionary)
    else:
        raise common.MFCException("Unrecognized input file format. Supported: .py, .json, .yaml, .yml. Please check the README and sample cases in the examples directory.")

    try:
        dictionary = json.loads(json_str)
    except Exception as exc:
        raise common.MFCException(f"Input file {filename} did not produce valid JSON. It should only print the case dictionary.\n\n{exc}\n")

    input_file = MFCInputFile(filename, dirpath, dictionary)
    input_file.validate_params(f"Input file {filename}")
    return input_file


load.CACHED_MFCInputFile = None
