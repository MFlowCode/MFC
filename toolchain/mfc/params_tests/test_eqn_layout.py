BASE = dict(
    num_fluids=2,
    num_vels=2,
    num_dims=2,
    nb_in=1,
    nmom_in=6,
    num_species=10,
    n=5,
    igr=False,
    bubbles_euler=False,
    qbmm=False,
    polytropic=True,
    adv_n=False,
    mhd=False,
    hypoelasticity=False,
    cyl_coord=False,
    surface_tension=False,
    cont_damage=False,
    hyper_cleaning=False,
    chemistry=False,
    six_eqn_alf_is_advected=False,
)


def layout(**kw):
    from mfc.params.eqn_layout import evaluate_layout

    return evaluate_layout({**BASE, "model_eqns": "5eq", **kw})


def test_five_equation_layout():
    out = layout()
    assert (out["cont"], out["mom"], out["E"], out["adv"]) == ((1, 2), (3, 4), 5, (6, 7))
    assert out["sys_size"] == 7 and out["alf"] == 1


def test_extensions_append_in_order():
    out = layout(bubbles_euler=True, polytropic=False, nb_in=3, hypoelasticity=True, chemistry=True)
    assert out["bub"] == (8, 19) and out["alf"] == 7
    assert out["stress"] == (20, 22) and out["species"] == (23, 32) and out["sys_size"] == 32


def test_grid_n_is_not_the_bubble_index():
    # eqn_idx%n (bubble number density) must not shadow the grid size n used by the MHD field count
    out = layout(bubbles_euler=True, adv_n=True, mhd=True, n=0)
    assert out["B"][1] - out["B"][0] + 1 == 2


def test_fortran_translation():
    from mfc.params.eqn_layout import _f, fortran_layout

    assert _f("adv.end if bubbles_euler else 1") == "merge(eqn_idx%adv%end, 1, bubbles_euler)"
    assert _f("num_dims*(num_dims + 1)//2") == "((num_dims * (num_dims + 1)) / 2)"
    assert "eqn_idx%species%end = sys_size + num_species" in fortran_layout()
