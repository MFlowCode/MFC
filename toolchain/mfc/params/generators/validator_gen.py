"""
Validator Code Generator.

Generates validation code from the constraint registry.
This allows maintaining validation rules in a single place while still
producing the check methods expected by CaseValidator.
"""

from typing import List
from ..schema import Stage, ConstraintType, ConstraintRule
from ..registry import ParamRegistry


def generate_validator_code(registry: ParamRegistry, stage: Stage) -> str:
    """
    Generate validation method code for a stage.

    Produces Python code that can be integrated into CaseValidator
    or used as a standalone validation function.

    Args:
        registry: The parameter registry with constraints
        stage: The stage to generate validation for

    Returns:
        Python source code as a string
    """
    constraints = registry.get_constraints_for_stage(stage)

    lines = [
        f'def validate_{stage.name.lower()}(self):',
        f'    """Auto-generated validation for {stage.name} stage."""',
    ]

    # Group constraints by check_method if specified
    by_method = {}
    for rule in constraints:
        method = rule.check_method or "general"
        if method not in by_method:
            by_method[method] = []
        by_method[method].append(rule)

    for method_name, rules in sorted(by_method.items()):
        if method_name != "general":
            lines.append(f'')
            lines.append(f'    # {method_name}')

        for rule in rules:
            code = _generate_constraint_check(rule)
            for line in code:
                lines.append(f'    {line}')

    if len(lines) == 2:
        lines.append('    pass')

    return '\n'.join(lines)


def _generate_constraint_check(rule: ConstraintRule) -> List[str]:
    """Generate code for a single constraint check."""
    lines = []

    if rule.constraint_type == ConstraintType.REQUIRED:
        param = rule.params[0]
        msg = rule.message or f"{param} is required"
        lines.append(f'self.require("{param}", "{msg}")')

    elif rule.constraint_type == ConstraintType.POSITIVE:
        param = rule.params[0]
        msg = rule.message or f"{param} must be positive"
        lines.append(f'val = self.get("{param}")')
        lines.append(f'if val is not None:')
        lines.append(f'    self.prohibit(val <= 0, "{msg}")')

    elif rule.constraint_type == ConstraintType.NON_NEGATIVE:
        param = rule.params[0]
        msg = rule.message or f"{param} must be non-negative"
        lines.append(f'val = self.get("{param}")')
        lines.append(f'if val is not None:')
        lines.append(f'    self.prohibit(val < 0, "{msg}")')

    elif rule.constraint_type == ConstraintType.RANGE:
        param = rule.params[0]
        lines.append(f'val = self.get("{param}")')
        lines.append(f'if val is not None:')

        if rule.min_value is not None:
            op = '<' if rule.min_exclusive else '<='
            inv_op = '<=' if rule.min_exclusive else '<'
            msg = rule.message or f"{param} must be > {rule.min_value}"
            lines.append(f'    self.prohibit(val {inv_op} {rule.min_value}, "{msg}")')

        if rule.max_value is not None:
            op = '>' if rule.max_exclusive else '>='
            inv_op = '>=' if rule.max_exclusive else '>'
            msg = rule.message or f"{param} must be < {rule.max_value}"
            lines.append(f'    self.prohibit(val {inv_op} {rule.max_value}, "{msg}")')

    elif rule.constraint_type == ConstraintType.CHOICES:
        param = rule.params[0]
        valid = rule.valid_values
        msg = rule.message or f"{param} must be one of {valid}"
        lines.append(f'val = self.get("{param}")')
        lines.append(f'if val is not None:')
        lines.append(f'    self.prohibit(val not in {valid}, "{msg}")')

    elif rule.constraint_type == ConstraintType.DEPENDENCY:
        lines.append(f'# Dependency: if {rule.if_param} {rule.if_condition}, '
                     f'then {rule.then_param} {rule.then_condition}')
        lines.append(f'if_val = self.get("{rule.if_param}")')
        lines.append(f'then_val = self.get("{rule.then_param}")')
        if_cond = _python_condition("if_val", rule.if_condition)
        then_cond = _python_condition("then_val", rule.then_condition)
        msg = rule.message
        lines.append(f'if {if_cond}:')
        lines.append(f'    self.prohibit(not ({then_cond}), "{msg}")')

    elif rule.constraint_type == ConstraintType.MUTEX:
        params = rule.mutex_params
        msg = rule.message or f"Only one of {params} can be enabled"
        lines.append(f'# Mutex: only one of {params}')
        lines.append(f'active = []')
        for p in params:
            lines.append(f'if self.get("{p}") in (True, "T", "t", 1):')
            lines.append(f'    active.append("{p}")')
        lines.append(f'self.prohibit(len(active) > 1, "{msg}")')

    elif rule.constraint_type == ConstraintType.CUSTOM:
        lines.append(f'# Custom constraint: {rule.rule_id}')
        lines.append(f'# Requires manual implementation')

    return lines


def _python_condition(var: str, condition: str) -> str:
    """Convert a condition string to Python expression."""
    if condition == "is not None":
        return f"{var} is not None"
    if condition == "is None":
        return f"{var} is None"
    if condition.startswith("> "):
        return f"{var} is not None and {var} > {condition[2:]}"
    if condition.startswith(">= "):
        return f"{var} is not None and {var} >= {condition[3:]}"
    if condition.startswith("< "):
        return f"{var} is not None and {var} < {condition[2:]}"
    if condition.startswith("<= "):
        return f"{var} is not None and {var} <= {condition[3:]}"
    if condition.startswith("== "):
        return f"{var} == {condition[3:]}"
    if condition.startswith("!= "):
        return f"{var} != {condition[3:]}"
    if condition == "is True" or condition == "== True":
        return f'{var} in (True, "T", "t", 1)'
    if condition == "is False" or condition == "== False":
        return f'{var} in (False, "F", "f", 0, None)'
    return f"bool({var})"


def generate_full_validator(registry: ParamRegistry) -> str:
    """
    Generate a complete validator module from the registry.

    This produces a standalone validator that can be used alongside
    or instead of the existing CaseValidator.

    Args:
        registry: The parameter registry with all constraints

    Returns:
        Python source code as a string
    """
    lines = [
        '"""',
        'Auto-generated Parameter Validator.',
        '',
        'Generated from mfc.params registry - Do not edit manually.',
        'To regenerate: python -m mfc.params.generators.validator_gen',
        '"""',
        '',
        'from typing import Dict, Any, List',
        '',
        '',
        'class GeneratedValidator:',
        '    """Validator with auto-generated constraint checks."""',
        '',
        '    def __init__(self, params: Dict[str, Any]):',
        '        self.params = params',
        '        self.errors: List[str] = []',
        '',
        '    def get(self, key: str, default=None):',
        '        """Get a parameter value."""',
        '        return self.params.get(key, default)',
        '',
        '    def require(self, key: str, message: str):',
        '        """Require a parameter to be set."""',
        '        if self.get(key) is None:',
        '            self.errors.append(message)',
        '',
        '    def prohibit(self, condition: bool, message: str):',
        '        """Add error if condition is True."""',
        '        if condition:',
        '            self.errors.append(message)',
        '',
    ]

    # Generate validation method for each stage
    for stage in [Stage.PRE_PROCESS, Stage.SIMULATION, Stage.POST_PROCESS]:
        lines.append('')
        lines.append(generate_validator_code(registry, stage))

    lines.append('')
    lines.append('    def validate_all(self) -> List[str]:')
    lines.append('        """Run all validation checks."""')
    lines.append('        self.validate_pre_process()')
    lines.append('        self.validate_simulation()')
    lines.append('        self.validate_post_process()')
    lines.append('        return self.errors')

    return '\n'.join(lines)


if __name__ == "__main__":
    from ..registry import REGISTRY
    from ..definitions import load_all_definitions
    load_all_definitions()

    print(generate_full_validator(REGISTRY))
