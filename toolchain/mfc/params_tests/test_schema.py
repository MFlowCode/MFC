"""Tests for ParamDef str_len field."""


def test_paramdef_str_len_default():
    from mfc.params.schema import ParamDef, ParamType

    p = ParamDef(name="foo", param_type=ParamType.STR)
    assert p.str_len == "name_len"


def test_paramdef_str_len_override():
    from mfc.params.schema import ParamDef, ParamType

    p = ParamDef(name="case_dir", param_type=ParamType.STR, str_len="path_len")
    assert p.str_len == "path_len"


def test_case_dir_has_path_len():
    import mfc.params.definitions  # noqa: F401
    from mfc.params.registry import REGISTRY

    p = REGISTRY.get_param_def("case_dir")
    assert p is not None
    assert p.str_len == "path_len"


def test_other_str_param_has_default_len():
    import mfc.params.definitions  # noqa: F401
    from mfc.params.registry import REGISTRY

    # cantera_file is another STR param that should use the default name_len
    p = REGISTRY.get_param_def("cantera_file")
    assert p is not None
    assert p.str_len == "name_len"
