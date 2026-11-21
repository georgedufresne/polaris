"""Optional primary-EEP rules (run/input_dsee.nml).

``input.nml`` leaves every rule at the Fortran behaviour; the golden test
already pins that.  These checks cover the switches themselves on the
synthetic track: the RGBBase primary, the He-ignition gate on the RGB tip, the
later-of-two tip rule, the very-low-mass fallback age, and the FSPS phases when
RGBBase is present.

Run directly (``python tests/test_dsee_options.py``) or under pytest.
"""

from __future__ import annotations

import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PY_ROOT = os.path.dirname(HERE)
RUN = os.path.join(PY_ROOT, 'run')
sys.path.insert(0, PY_ROOT)
sys.path.insert(0, HERE)

import make_eep  # noqa: E402
from make_synthetic_track import (  # noqa: E402
    write_history_file,
    write_history_columns_list,
)
from pyeep.eep import RGBTip, TAMS, primary_eep, secondary_eep  # noqa: E402
from pyeep.phases import RGB, RED_HB_CLUMP, set_track_phase  # noqa: E402
from pyeep.support import (  # noqa: E402
    EEPState,
    Track,
    read_history_file,
    set_eep_interval,
    setup_columns,
)


def _state(nml: str, eep: str) -> EEPState:
    state = EEPState()
    make_eep.apply_eep_controls(make_eep.read_namelist_group(nml, 'eep_controls'), state)
    set_eep_interval(eep, state)
    return state


def _dsee_state() -> EEPState:
    return _state(os.path.join(RUN, 'input_dsee.nml'), os.path.join(RUN, 'input_dsee.eep'))


def _default_state() -> EEPState:
    return _state(os.path.join(RUN, 'input.nml'), os.path.join(RUN, 'input.eep'))


def _track(state: EEPState, initial_mass: float = 1.0, mutate=None) -> Track:
    """Synthetic track read through the normal history reader, optionally edited."""
    with tempfile.TemporaryDirectory() as root:
        write_history_file(os.path.join(root, 'x.track'), initial_mass=initial_mass)
        cols = os.path.join(root, 'cols.list')
        write_history_columns_list(cols)
        setup_columns(cols, state)
        t = Track()
        t.filename = 'x.track'
        t.ncol = state.ncol
        t.neep = state.primary
        read_history_file(t, state, root)
    if mutate is not None:
        mutate(t, state)
    return t


def _with_rgb(t: Track, state: EEPState) -> None:
    """Keep central He at its post-MS value up to the luminosity peak, so the
    track has a giant branch (in the raw synthetic track He burning starts
    right at TAMS and the RGB-tip search window is a few rows long)."""
    Yc = t.tr[state.i_Yc]
    peak = int(np.argmax(t.tr[state.i_logL]))
    tams = int(np.argmax(t.tr[state.i_Xc] < 1e-12))
    Yc[tams:peak + 1] = Yc[tams]
    Yc[peak + 1:] = np.minimum(Yc[peak + 1:], Yc[tams] - 0.02)


def _with_knee(t: Track, state: EEPState) -> None:
    """Make TAMS -> tip a subgiant branch (cooling at nearly constant L) then a
    giant branch (rising L at nearly constant Teff), with the turn at row KNEE."""
    _with_rgb(t, state)
    tams = int(np.argmax(t.tr[state.i_Xc] < 1e-12)) + 1
    tip = int(np.argmax(t.tr[state.i_logL])) + 1
    knee = tams + (tip - tams) // 3
    T, L = t.tr[state.i_logTe], t.tr[state.i_logL]
    a, k, b = tams - 1, knee - 1, tip - 1
    T[a:k + 1] = np.linspace(T[a], 3.62, k - a + 1)
    L[a:k + 1] = np.linspace(L[a], L[a] + 0.1, k - a + 1)
    T[k:b + 1] = np.linspace(3.62, 3.55, b - k + 1)
    L[k:b + 1] = np.linspace(L[a] + 0.1, L[b], b - k + 1)
    t._knee = knee


def test_dsee_namelist_switches_everything_on():
    d, s = _default_state(), _dsee_state()
    assert (s.IAMS_Xc, s.RGBTip_require_He_ignition, s.RGBTip_take_later,
            s.add_RGBBase, s.merger_fallback) == (0.30, True, True, True, False)
    assert s.logL_scale == 2.0 and s.VLM_fallback_min_age == 2.09e10
    assert (d.IAMS_Xc, d.RGBTip_require_He_ignition, d.RGBTip_take_later,
            d.add_RGBBase, d.merger_fallback) == (0.35, False, False, False, True)
    assert d.VLM_fallback_min_age == 1.5e10
    # .eep rows of the primaries: TAMS 454 and RGBTip 605 in both, RGBBase 491
    rows = lambda st: np.cumsum([1] + [n + 1 for n in st.eep_interval])
    assert list(rows(d)[:5]) == [1, 202, 353, 454, 605]
    assert list(rows(s)[:6]) == [1, 202, 353, 454, 491, 605]


def test_v11_namelist():
    """input_dsee_v1.1.nml: the 1.1 rules on the standard metric and EEP intervals (no RGBBase)."""
    s = _state(os.path.join(RUN, 'input_dsee_v1.1.nml'), os.path.join(RUN, 'input.eep'))
    assert (s.IAMS_Xc, s.RGBTip_require_He_ignition, s.RGBTip_take_later,
            s.add_RGBBase, s.merger_fallback) == (0.30, True, True, False, False)
    assert s.logL_scale == 0.125 and s.VLM_fallback_min_age == 1.5e10 and s.primary == 9
    nml = make_eep.read_namelist_group(os.path.join(RUN, 'input_dsee_v1.1.nml'), 'eep_controls')
    assert str(nml['eep_input_file']).strip("'\" ") == 'input.eep'


def test_rgbbase_sits_at_the_knee():
    state = _dsee_state()
    t = _track(state, mutate=_with_knee)
    primary_eep(t, state)
    tams, base, tip = int(t.eep[3]), int(t.eep[4]), int(t.eep[5])
    assert tams < base < tip, t.eep
    assert abs(base - t._knee) <= 2, (base, t._knee)
    # the tip is the same row the rule finds without a base: the luminosity peak
    assert tip == RGBTip(t, state, tams + 1) == int(np.argmax(t.tr[state.i_logL])) + 1


def test_rgb_tip_needs_he_ignition():
    def no_ignition(t, state):
        _with_rgb(t, state)
        t.tr[state.i_logLHe][:] = -30.0
    state = _dsee_state()
    t = _track(state, mutate=no_ignition)
    primary_eep(t, state)
    assert t.eep[3] > 0 and t.eep[4] == 0 and t.eep[5] == 0, t.eep  # TAMS, no base, no tip
    state.RGBTip_require_He_ignition = False
    t = _track(state, mutate=no_ignition)
    assert RGBTip(t, state, 2) > 0


def test_rgb_tip_later_rule():
    def cool_after_peak(t, state):
        i = int(np.argmax(t.tr[state.i_logL]))
        _with_rgb(t, state)
        t.tr[state.i_Yc][i + 1:i + 4] = t.tr[state.i_Yc][i]                # still on the RGB 3 rows later
        t.tr[state.i_logTe][i + 3] = t.tr[state.i_logTe][i] - 0.01   # coolest point 3 rows after max L
    state = _default_state()
    t = _track(state, mutate=cool_after_peak)
    tams = TAMS(t, state, 1e-12, 1)
    early = RGBTip(t, state, tams + 1)
    state.RGBTip_take_later = True
    late = RGBTip(t, state, tams + 1)
    assert late == early + 3, (early, late)


def test_vlm_fallback_needs_min_age():
    def never_leaves_ms(t, state):
        t.tr[state.i_Xc][:] = np.maximum(t.tr[state.i_Xc], 0.5)
    state = _default_state()
    state.very_low_mass_limit = 0.7
    age_end = None
    for min_age, expect in ((1.5e10, 0), (1.0e10, 'last')):
        state.VLM_fallback_min_age = min_age
        t = _track(state, initial_mass=0.5, mutate=never_leaves_ms)
        age_end = t.tr[state.i_age][-1]
        got = TAMS(t, state, 1e-12, 2)
        if expect == 0:
            assert got == 0 and not t.vlm_fallback
        else:
            assert got == t.ntrack and t.vlm_fallback
    assert 1.0e10 <= age_end < 1.5e10


def test_phases_with_rgbbase():
    state = _dsee_state()
    t = _track(state, mutate=_with_knee)
    primary_eep(t, state)
    s = secondary_eep(t, state)
    set_track_phase(s, state)
    assert np.all(s.phase[453:604] == RGB)                # EEP 454-604: TAMS -> RGBBase -> RGBTip
    assert s.phase[604] == RED_HB_CLUMP                   # the tip row, as in the Fortran
    # same phase on every .eep row as without RGBBase
    d = _default_state()
    t0 = _track(d, mutate=_with_knee)
    primary_eep(t0, d)
    s0 = secondary_eep(t0, d)
    set_track_phase(s0, d)
    n = min(s.ntrack, s0.ntrack)
    assert n > 605 and np.array_equal(s.phase[:n], s0.phase[:n])


def test_pipeline_with_dsee_rules():
    """make_eep with input_dsee.nml/.eep, then make_iso: the RGBBase primary is
    found on every mass and the isochrones still recover the analytic masses."""
    import make_iso
    import test_pipeline as tp
    from pyeep.iso import read_eep
    ctrl = open(os.path.join(RUN, 'input_dsee.nml')).read()
    ctrl = ctrl[ctrl.index('&eep_controls'):]
    ctrl = ctrl[:ctrl.index('/ !end of eep_controls')] + '/\n'
    ctrl = ctrl.replace("do_phases = .false.", "do_phases = .true.")
    with tempfile.TemporaryDirectory() as root:
        os.makedirs(os.path.join(root, 'hist'))
        names = []
        for m in tp.MASSES:
            names.append(f"m{int(round(m * 100)):04d}")
            tp._write_history(os.path.join(root, 'hist', names[-1]), m)
        write_history_columns_list(os.path.join(root, 'cols.list'))
        with open(os.path.join(RUN, 'input_dsee.eep')) as fh:
            eep_text = fh.read()
        for fname, text in (('input_dsee.eep', eep_text), ('input.format', tp.INPUT_FORMAT),
                            ('input.nml', ctrl + '\n&iso_controls\n   do_linear_interpolation = .true.\n'
                                          '   log_age_delta = 1.0\n/\n')):
            with open(os.path.join(root, fname), 'w') as fh:
                fh.write(text)
        with open(os.path.join(root, 'input.example'), 'w') as fh:
            fh.write(f"#v\npipe\n#Y Z FeH aFe v\n 0.2700 1.40000e-02 0.00 0.00 0.0\n"
                     f"#dirs\n{root}/hist\n{root}/eeps\n{root}/isos\n#cols\n{root}/cols.list\n"
                     f"#tracks\n{len(names)}\n" + ''.join(n + '\n' for n in names)
                     + f"#isochrones\npipe.iso\nlist\nlinear\n{len(tp.AGES)}\n"
                     + ''.join(f"{a:.6e}\n" for a in tp.AGES))
        cwd = os.getcwd()
        try:
            os.chdir(root)
            assert make_eep.run('input.example') == 0
            isos = make_iso.make_isochrones('input.example')
        finally:
            os.chdir(cwd)
        tracks = [read_eep(os.path.join(root, 'eeps', n + '.eep')) for n in names]
    ref = tracks[0]
    assert ref.neep >= 6 and ref.eep[3] < ref.eep[4] < ref.eep[5], ref.eep   # TAMS < RGBBase < RGBTip
    for t in tracks:
        assert t.eep.tolist() == ref.eep.tolist()
    age_ref = ref.tr[ref.col_names.index('star_age')]
    for iso, age in zip(isos, tp.AGES):
        mi = iso.column('initial_mass')
        assert iso.neep > 20 and np.all(np.diff(mi) >= -1e-12)
        m_true = tp.MASSES[0] * (age / age_ref[iso.eep - 1]) ** (-1 / tp.SLOPE)
        assert np.max(np.abs(mi - m_true)) < 5e-3


def _ladder_state() -> EEPState:
    return _state(os.path.join(RUN, 'input_dsee_ladder.nml'),
                  os.path.join(RUN, 'input_dsee_ladder.eep'))


def test_ladder_namelist_and_intervals():
    """input_dsee_ladder.nml/.eep: 19 primaries, TAMS on 454 and the last on 605."""
    s = _ladder_state()
    assert s.Xc_HeL_ladder and s.primary == 19 and len(s.eep_interval) == 18
    assert not _dsee_state().Xc_HeL_ladder and not _default_state().Xc_HeL_ladder
    rows = np.cumsum([1] + [n + 1 for n in s.eep_interval])
    assert list(rows) == [1, 202, 213, 224, 243, 262, 281, 300, 319, 338, 363,
                          394, 429, 454, 479, 492, 517, 554, 605], list(rows)
    nml = make_eep.read_namelist_group(os.path.join(RUN, 'input_dsee_ladder.nml'),
                                       'eep_controls')
    assert str(nml['eep_input_file']).strip("\'\" ") == 'input_dsee_ladder.eep'


def test_ladder_primaries():
    """The ladder places 19 strictly increasing primaries on the synthetic track."""
    state = _ladder_state()
    t = _track(state, mutate=_with_rgb)
    primary_eep(t, state)
    eep = t.eep.tolist()
    assert len(eep) == 19 and all(e > 0 for e in eep), eep
    assert all(b > a for a, b in zip(eep, eep[1:])), eep

    # rungs 3-13 sit where Xc has fallen to the requested fraction of Xc(ZAMS),
    # and the row before each is still above it
    Xc = t.tr[state.i_Xc]
    Xc_zams = Xc[eep[1] - 1]
    for k, frac in enumerate(state.Xc_ladder):
        i = eep[2 + k] - 1
        assert Xc[i] <= frac * Xc_zams < Xc[i - 1], (frac, Xc[i])

    # TAMS is the usual Xc = 1e-12 anchor, with no IAMS before it
    assert eep[13] == TAMS(t, state, 1e-12, eep[12] + 1)

    # rungs 15-19: helium luminosity crossing each decade, all after the TAMS
    ratio = t.tr[state.i_logLHe] - t.tr[state.i_logL]
    for k, level in enumerate(state.HeL_ladder):
        i = eep[14 + k] - 1
        assert i >= eep[13] and ratio[i] >= level, (level, ratio[i])


def test_ladder_stops_where_the_track_does():
    """A track whose core never burns down past a rung keeps the primaries below
    it and leaves the rest at 0, rather than piling them on the last row."""
    state = _ladder_state()
    floor = [0.0]

    def _young(t, st):
        Xc = t.tr[st.i_Xc]
        floor[0] = 0.55 * Xc[0]
        Xc[:] = np.maximum(Xc, floor[0])

    t = _track(state, mutate=_young)
    primary_eep(t, state)
    eep = t.eep.tolist()
    Xc_zams = t.tr[state.i_Xc][eep[1] - 1]
    reached = sum(1 for f in state.Xc_ladder if f * Xc_zams >= floor[0])
    assert 0 < reached < len(state.Xc_ladder), reached
    placed = 2 + reached
    assert all(e > 0 for e in eep[:placed]), eep
    assert all(b > a for a, b in zip(eep[:placed], eep[1:placed])), eep
    assert all(e == 0 for e in eep[placed:]), eep


if __name__ == '__main__':
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            fn()
            print(name, 'OK')
