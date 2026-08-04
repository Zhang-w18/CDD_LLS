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
