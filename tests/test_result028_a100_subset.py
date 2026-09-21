from tools import plot_result028_a100_subset as subset


def test_fft_sample_delay_conversion_preserves_physical_delay() -> None:
    assert subset.fft_sample_delays([0.0, 9.0, 63.0]) == [0.0, 64.0, 448.0]


def test_subset_rows_keeps_selected_candidates_through_16db() -> None:
    rows = [
        {"candidate_id": "A100_B0_QC", "snr_db": "16.0"},
        {"candidate_id": "A100_B0_QC", "snr_db": "16.25"},
        {"candidate_id": "A100_AP_TEPS_T1", "snr_db": "15.0"},
        {"candidate_id": "A100_MEFF_T2_06", "snr_db": "14.0"},
    ]
    kept = subset.subset_rows(rows)
    assert [(row["candidate_id"], row["snr_db"]) for row in kept] == [
        ("A100_B0_QC", "16.0"),
        ("A100_MEFF_T2_06", "14.0"),
    ]


def test_json_dash_style_is_converted_for_matplotlib() -> None:
    assert subset._matplotlib_linestyle([0, [1, 1]]) == (0, (1, 1))


def test_subset_rows_can_extend_small_delay_curve_through_21db() -> None:
    rows = [
        {"candidate_id": "A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD", "snr_db": "21.0"},
        {"candidate_id": "A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD", "snr_db": "21.25"},
    ]
    kept = subset.subset_rows(rows, subset.SMALL_DELAY_SELECTED, snr_max_db=21.0)
    assert [row["snr_db"] for row in kept] == ["21.0"]


def test_ideal_subset_can_respect_statistical_plot_cutoff() -> None:
    rows = [
        {
            "candidate_id": "A100_S0_SIDON",
            "snr_db": "15.5",
            "plot_included": "True",
        },
        {
            "candidate_id": "A100_S0_SIDON",
            "snr_db": "15.75",
            "plot_included": "False",
        },
    ]
    kept = subset.subset_rows(rows, require_plot_included=True)
    assert [row["snr_db"] for row in kept] == ["15.5"]


def test_section37_uses_requested_receiver_comparison_and_renames_delay_set_2() -> None:
    assert subset.SECTION37_SELECTED == (
        ("large-delay CDD, transparent", "A100_AP_RMS_T1_TRANSPARENT_CDD"),
        ("large-delay CDD, non-transparent", "A100_AP_RMS_T1"),
        ("small-delay CDD, transparent", "A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD"),
        ("small-delay CDD, non-transparent", "A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD"),
        ("precoder cycling 6 RB, transparent", "A100_PRG_DFT8_6RB"),
    )
    assert subset.SECTION37_X_MIN_DB == 13.75
    assert subset.SECTION37_X_MAX_DB == 20.25
    assert subset.SECTION37_SNR_MAX_BY_CANDIDATE == {
        "A100_PRG_DFT8_6RB": 17.0,
        "A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD": 18.5,
    }
