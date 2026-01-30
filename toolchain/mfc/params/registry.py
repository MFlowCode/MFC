"""
Parameter Registry.

Central coordinator for all MFC parameter definitions and constraints.
This is the single source of truth that generators use to produce
case_dicts.py schemas and validator code.
"""

from typing import Dict, List, Optional, Set, Any
from collections import defaultdict

from .schema import (
    ParamDef, ConstraintRule, PatternDef,
    ParamType, Stage, ConstraintType
)


class ParamRegistry:
    """
    Central registry for MFC parameters and constraints.

    Provides:
    - Registration of parameters and patterns
    - Registration of validation constraints
    - Query methods for parameters by stage/category
    - Query methods for constraints by stage/type
    - Export methods for code generation
    """

    def __init__(self):
        self._params: Dict[str, ParamDef] = {}
        self._patterns: List[PatternDef] = []
        self._constraints: Dict[str, ConstraintRule] = {}
        self._expanded = False

        # Index for faster lookups
        self._params_by_stage: Dict[Stage, Set[str]] = defaultdict(set)
        self._params_by_category: Dict[str, Set[str]] = defaultdict(set)
        self._constraints_by_stage: Dict[Stage, List[str]] = defaultdict(list)
        self._constraints_by_type: Dict[ConstraintType, List[str]] = defaultdict(list)

    def register(self, param: ParamDef, merge_stages: bool = True) -> None:
        """
        Register a parameter definition.

        Args:
            param: ParamDef instance to register
            merge_stages: If True and param already exists, merge stages instead of error

        Raises:
            ValueError: If parameter exists with different type and merge_stages=False
        """
        if param.name in self._params:
            existing = self._params[param.name]
            if merge_stages:
                # Merge stages - parameter can exist in multiple stages
                if existing.param_type != param.param_type:
                    raise ValueError(
                        f"Parameter '{param.name}' already registered with type "
                        f"{existing.param_type}, cannot register with type {param.param_type}"
                    )
                # Add new stages to existing parameter
                existing.stages.update(param.stages)
                for stage in param.stages:
                    self._params_by_stage[stage].add(param.name)
                return
            raise ValueError(f"Parameter '{param.name}' already registered")

        self._params[param.name] = param

        # Update indices
        for stage in param.stages:
            self._params_by_stage[stage].add(param.name)
        if param.category:
            self._params_by_category[param.category].add(param.name)

    def register_pattern(self, pattern: PatternDef) -> None:
        """
        Register a parameterized pattern for bulk parameter generation.

        Patterns are expanded into concrete parameters when expand_patterns() is called
        or when querying parameters.

        Args:
            pattern: PatternDef instance to register
        """
        self._patterns.append(pattern)
        self._expanded = False

    def register_constraint(self, rule: ConstraintRule) -> None:
        """
        Register a validation constraint.

        Args:
            rule: ConstraintRule instance to register

        Raises:
            ValueError: If constraint with same rule_id already registered
        """
        if rule.rule_id in self._constraints:
            raise ValueError(f"Constraint '{rule.rule_id}' already registered")

        self._constraints[rule.rule_id] = rule

        # Update indices
        for stage in rule.stages:
            self._constraints_by_stage[stage].append(rule.rule_id)
        self._constraints_by_type[rule.constraint_type].append(rule.rule_id)

    def expand_patterns(self, index_limits: Optional[Dict[str, int]] = None) -> None:
        """
        Expand all registered patterns into concrete parameters.

        Args:
            index_limits: Dict mapping index name to max value.
                          If None, uses default limits.
        """
        if index_limits is None:
            # Default limits matching MFC's case_dicts.py
            index_limits = {
                "num_patches": 10,
                "num_fluids": 10,
                "num_probes": 10,
                "num_source": 5,
                "num_ibs": 5,
                "num_bc_patches": 10,
                "Nb": 5,  # bubble bins
            }

        for pattern in self._patterns:
            for param in pattern.expand(index_limits):
                if param.name not in self._params:
                    self._params[param.name] = param
                    for stage in param.stages:
                        self._params_by_stage[stage].add(param.name)
                    if param.category:
                        self._params_by_category[param.category].add(param.name)

        self._expanded = True

    def get(self, name: str) -> Optional[ParamDef]:
        """Get a parameter by name."""
        return self._params.get(name)

    def get_constraint(self, rule_id: str) -> Optional[ConstraintRule]:
        """Get a constraint by rule_id."""
        return self._constraints.get(rule_id)

    def get_params_by_stage(self, stage: Stage) -> Dict[str, ParamDef]:
        """
        Get all parameters valid for a given stage.

        Args:
            stage: The stage to query (COMMON, PRE_PROCESS, etc.)

        Returns:
            Dict mapping parameter names to ParamDef instances
        """
        result = {}
        # Include COMMON parameters for all stages
        for name in self._params_by_stage.get(Stage.COMMON, set()):
            result[name] = self._params[name]
        if stage != Stage.COMMON:
            for name in self._params_by_stage.get(stage, set()):
                result[name] = self._params[name]
        return result

    def get_params_by_category(self, category: str) -> Dict[str, ParamDef]:
        """Get all parameters in a given category."""
        return {
            name: self._params[name]
            for name in self._params_by_category.get(category, set())
        }

    def get_constraints_for_stage(self, stage: Stage) -> List[ConstraintRule]:
        """
        Get all constraints that apply to a given stage.

        Args:
            stage: The stage to query

        Returns:
            List of ConstraintRule instances, sorted by priority
        """
        result = []
        seen = set()

        # Include COMMON constraints for all stages
        for rule_id in self._constraints_by_stage.get(Stage.COMMON, []):
            if rule_id not in seen:
                result.append(self._constraints[rule_id])
                seen.add(rule_id)

        if stage != Stage.COMMON:
            for rule_id in self._constraints_by_stage.get(stage, []):
                if rule_id not in seen:
                    result.append(self._constraints[rule_id])
                    seen.add(rule_id)

        return sorted(result, key=lambda r: r.priority)

    def get_constraints_by_type(self, ctype: ConstraintType) -> List[ConstraintRule]:
        """Get all constraints of a given type."""
        return [
            self._constraints[rule_id]
            for rule_id in self._constraints_by_type.get(ctype, [])
        ]

    @property
    def all_params(self) -> Dict[str, ParamDef]:
        """Get all registered parameters."""
        return dict(self._params)

    @property
    def all_constraints(self) -> Dict[str, ConstraintRule]:
        """Get all registered constraints."""
        return dict(self._constraints)

    @property
    def param_count(self) -> int:
        """Total number of registered parameters."""
        return len(self._params)

    @property
    def constraint_count(self) -> int:
        """Total number of registered constraints."""
        return len(self._constraints)

    def to_case_dicts_schema(self, stage: Stage) -> Dict[str, str]:
        """
        Generate case_dicts.py compatible schema for a stage.

        Returns a dict mapping parameter names to type tags.
        """
        params = self.get_params_by_stage(stage)
        return {name: param.type_tag for name, param in params.items()}

    def validate(self, params: Dict[str, Any], stage: Stage) -> List[str]:
        """
        Validate parameters against registered constraints.

        Args:
            params: Dict of parameter values to validate
            stage: The stage being validated

        Returns:
            List of error messages (empty if all valid)
        """
        errors = []
        constraints = self.get_constraints_for_stage(stage)

        for rule in constraints:
            error = self._check_constraint(rule, params)
            if error:
                errors.append(error)

        return errors

    def _check_constraint(self, rule: ConstraintRule, params: Dict[str, Any]) -> Optional[str]:
        """Check a single constraint against parameters."""

        if rule.constraint_type == ConstraintType.REQUIRED:
            param_name = rule.params[0]
            if params.get(param_name) is None:
                return rule.message or f"{param_name} is required"

        elif rule.constraint_type == ConstraintType.POSITIVE:
            param_name = rule.params[0]
            value = params.get(param_name)
            if value is not None and isinstance(value, (int, float)):
                if value <= 0:
                    return rule.message or f"{param_name} must be positive (got {value})"

        elif rule.constraint_type == ConstraintType.NON_NEGATIVE:
            param_name = rule.params[0]
            value = params.get(param_name)
            if value is not None and isinstance(value, (int, float)):
                if value < 0:
                    return rule.message or f"{param_name} must be non-negative (got {value})"

        elif rule.constraint_type == ConstraintType.RANGE:
            param_name = rule.params[0]
            value = params.get(param_name)
            if value is not None and isinstance(value, (int, float)):
                if rule.min_value is not None:
                    if rule.min_exclusive:
                        if value <= rule.min_value:
                            return rule.message or f"{param_name} must be > {rule.min_value}"
                    else:
                        if value < rule.min_value:
                            return rule.message or f"{param_name} must be >= {rule.min_value}"
                if rule.max_value is not None:
                    if rule.max_exclusive:
                        if value >= rule.max_value:
                            return rule.message or f"{param_name} must be < {rule.max_value}"
                    else:
                        if value > rule.max_value:
                            return rule.message or f"{param_name} must be <= {rule.max_value}"

        elif rule.constraint_type == ConstraintType.CHOICES:
            param_name = rule.params[0]
            value = params.get(param_name)
            if value is not None and value not in rule.valid_values:
                return rule.message or f"{param_name} must be one of {rule.valid_values}"

        elif rule.constraint_type == ConstraintType.DEPENDENCY:
            if_value = params.get(rule.if_param)
            then_value = params.get(rule.then_param)

            # Evaluate if_condition
            if_met = self._eval_condition(if_value, rule.if_condition)
            if if_met:
                # Check then_condition
                then_met = self._eval_condition(then_value, rule.then_condition)
                if not then_met:
                    return rule.message or (
                        f"When {rule.if_param} {rule.if_condition}, "
                        f"{rule.then_param} must {rule.then_condition}"
                    )

        elif rule.constraint_type == ConstraintType.MUTEX:
            active = [
                p for p in rule.mutex_params
                if params.get(p) in (True, 'T', 't', 1)
            ]
            if len(active) > 1:
                return rule.message or f"Only one of {rule.mutex_params} can be enabled"

        elif rule.constraint_type == ConstraintType.CUSTOM:
            if rule.predicate and not rule.predicate(params):
                return rule.message or f"Custom constraint {rule.rule_id} failed"

        return None

    def _eval_condition(self, value: Any, condition: str) -> bool:
        """Evaluate a simple condition string against a value."""
        if condition == "is not None":
            return value is not None
        if condition == "is None":
            return value is None
        if condition.startswith("> "):
            threshold = float(condition[2:])
            return value is not None and value > threshold
        if condition.startswith(">= "):
            threshold = float(condition[3:])
            return value is not None and value >= threshold
        if condition.startswith("< "):
            threshold = float(condition[2:])
            return value is not None and value < threshold
        if condition.startswith("<= "):
            threshold = float(condition[3:])
            return value is not None and value <= threshold
        if condition.startswith("== "):
            expected = condition[3:]
            if expected.startswith("'") and expected.endswith("'"):
                expected = expected[1:-1]
            elif expected.isdigit():
                expected = int(expected)
            return value == expected
        if condition.startswith("!= "):
            expected = condition[3:]
            if expected.startswith("'") and expected.endswith("'"):
                expected = expected[1:-1]
            elif expected.isdigit():
                expected = int(expected)
            return value != expected
        if condition == "is True" or condition == "== True":
            return value in (True, 'T', 't', 1)
        if condition == "is False" or condition == "== False":
            return value in (False, 'F', 'f', 0, None)

        # Default: treat as truth check
        return bool(value)

    def summary(self) -> str:
        """Return a summary of registered parameters and constraints."""
        lines = [
            "Parameter Registry Summary",
            "=" * 40,
            f"Total parameters: {self.param_count}",
            f"Total constraints: {self.constraint_count}",
            "",
            "Parameters by stage:",
        ]

        for stage in Stage:
            count = len(self._params_by_stage.get(stage, set()))
            lines.append(f"  {stage.name}: {count}")

        lines.append("")
        lines.append("Constraints by type:")

        for ctype in ConstraintType:
            count = len(self._constraints_by_type.get(ctype, []))
            if count > 0:
                lines.append(f"  {ctype.name}: {count}")

        return "\n".join(lines)


# Global registry instance
REGISTRY = ParamRegistry()
