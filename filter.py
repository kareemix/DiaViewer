def filter_diann(
    df,
    empirical_lib=False,
    peptidoform_mode=False,
    plexdia=False,
    PGMaxLFQ=False,
    QQ=False,
    avg_quality_filter=False,
    filter_peak_width=False,
    min_points_across_peak=8,
    duty_cycle=1.5  # in seconds
):
    def report(name, before, after):
        removed = before - after
        kept_pct = after / before * 100 if before > 0 else 0
        removed_pct = removed / before * 100 if before > 0 else 0
        print(f"{name}: {after} rows kept ({kept_pct:.2f}%), {removed} removed ({removed_pct:.2f}%)")

    original_count = len(df)
    print(f"Starting with {original_count} rows")

    # Step 1
    df1 = df[(df["Q.Value"] < 0.01) & (df["Lib.Q.Value"] < 0.01)]
    report("Q.Value & Lib.Q.Value < 0.01", original_count, len(df1))

    # Step 2
    if not empirical_lib:
        df2 = df1[
            (df1["Global.Q.Value"] < 0.01) &
            (df1["Global.PG.Q.Value"] < 0.01)
        ]
        report("Global Q filters", len(df1), len(df2))
    else:
        df2 = df1
        print("Skipped Global.Q filters (empirical lib)")

    # Step 3
    df3 = df2[
        (df2["Lib.PG.Q.Value"] < 0.01) &
        (df2["PG.Q.Value"] <= 0.05)
    ]
    report("PG filters", len(df2), len(df3))

    # Step 4
    if peptidoform_mode:
        df4 = df3[
            (df3["Lib.Peptidoform.Q.Value"] < 0.01) &
            (df3["Global.Peptidoform.Q.Value"] < 0.01)
        ]
        report("Peptidoform filters", len(df3), len(df4))
    else:
        df4 = df3
        print("Skipped Peptidoform filters")

    # Step 5
    if plexdia:
        df5 = df4[df4["Channel.Q.Value"] <= 0.5]
        report("Channel.Q filter", len(df4), len(df5))
    else:
        df5 = df4
        print("Skipped Channel.Q filter")

    # Step 6
    if QQ:
        df6 = df5[df5["Quantity.Quality"] >= 0.5]
        report("Quantity.Quality", len(df5), len(df6))
    else:
        df6 = df5
        print("Skipped Quantity.Quality filter")

     # Step 7
    if PGMaxLFQ:
        df7 = df6[df6["PG.MaxLFQ.Quality"] >= 0.7]
        report("PG.MaxLFQ.Quality", len(df6), len(df7))
    else:
        df7 = df6
        print("Skipped PG.MaxLFQ.Quality filter")

    # Step 8
    if avg_quality_filter:
        avg_q = df7.groupby("Precursor.Id")[["Quantity.Quality", "PG.MaxLFQ.Quality"]].mean()
        good_precursors = avg_q[
            (avg_q["Quantity.Quality"] >= 0.5) &
            (avg_q["PG.MaxLFQ.Quality"] >= 0.7)
        ].index
        df8 = df7[df7["Precursor.Id"].isin(good_precursors)]
        report("Avg Precursor.Id quality filter", len(df7), len(df8))
    else:
        df8 = df7
        print("Skipped avg Precursor.Id quality filter")

    # Step 9: Peak width filter
    if filter_peak_width:
        df8 = df8.copy()  # :bulb: This avoids SettingWithCopyWarning
        df8["PeakWidth_sec"] = (df8["RT.Stop"] - df8["RT.Start"]) * 60
        df8["PointsAcrossPeak"] = df8["PeakWidth_sec"] / duty_cycle
        df9 = df8[df8["PointsAcrossPeak"] >= min_points_across_peak]
        report(f"Peak Width ≥ {min_points_across_peak} points", len(df8), len(df9))
    else:
        df9 = df8
        print("Skipped Peak Width filter")

    final_count = len(df9)
    total_pct = final_count / original_count * 100 if original_count > 0 else 0
    print(f"Final count: {final_count} rows ({total_pct:.2f}% kept, {original_count - final_count} removed)")

    return df9


# df_nofilter = filter_diann(
#     df,
#     empirical_lib=True,
#     peptidoform_mode=True,
#     plexdia=False,
#     PGMaxLFQ=False,
#     QQ=False,
#     avg_quality_filter=False,
#     filter_peak_width=False,
#     min_points_across_peak=8,
#     duty_cycle=1.5  # in seconds
# )
