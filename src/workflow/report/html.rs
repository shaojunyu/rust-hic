//! Render an offline, self-contained report without external assets.
use super::{Report, Resolution};
use std::fmt::Write;

fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            _ => out.push(c),
        }
    }
    out
}

fn number(n: u64) -> String {
    let raw = n.to_string();
    let mut result = String::with_capacity(raw.len() + raw.len() / 3);
    for (i, c) in raw.chars().enumerate() {
        if i != 0 && (raw.len() - i).is_multiple_of(3) {
            result.push(',');
        }
        result.push(c);
    }
    result
}

fn percentage(n: u64, d: u64) -> String {
    if d == 0 {
        "—".into()
    } else {
        format!("{:.1}%", 100.0 * n as f64 / d as f64)
    }
}

fn duration(seconds: f64) -> String {
    if seconds < 60.0 {
        return format!("{seconds:.1} s");
    }
    let s = seconds.round() as u64;
    if s < 3600 {
        format!("{} min {:02} s", s / 60, s % 60)
    } else {
        format!("{} h {:02} min {:02} s", s / 3600, (s / 60) % 60, s % 60)
    }
}

fn resolution(n: u64) -> String {
    if n >= 1_000_000 && n.is_multiple_of(1_000_000) {
        format!("{} Mb", n / 1_000_000)
    } else if n >= 1_000_000 && n.is_multiple_of(100_000) {
        format!("{:.1} Mb", n as f64 / 1_000_000.0)
    } else if n >= 1_000 && n.is_multiple_of(1_000) {
        format!("{} kb", number(n / 1_000))
    } else {
        format!("{} bp", number(n))
    }
}

fn balance_state(r: &Resolution) -> (&'static str, &'static str) {
    match r.converged {
        Some(true) => ("Converged", "good"),
        Some(false) => ("Not converged", "warn"),
        None => ("No usable bins", "muted"),
    }
}

fn card(out: &mut String, label: &str, value: &str, detail: &str) {
    let _ = write!(out, "<article class=\"card\"><div class=\"eyebrow\">{}</div><div class=\"metric\">{}</div><p>{}</p></article>", escape(label), escape(value), escape(detail));
}

fn row(out: &mut String, label: &str, n: u64, denominator: Option<u64>, detail: &str) {
    let share = denominator.map_or("—".into(), |d| percentage(n, d));
    let _ = write!(out, "<tr><th scope=\"row\">{}</th><td class=\"numeric\">{}</td><td class=\"numeric\">{}</td><td>{}</td></tr>", escape(label), number(n), share, escape(detail));
}

pub(super) fn render(r: &Report) -> String {
    let mut out = String::with_capacity(20_000);
    out.push_str("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"><meta name=\"color-scheme\" content=\"light\"><title>rust-hic · Run report</title><style>");
    out.push_str(include_str!("style.css"));
    out.push_str("</style></head><body><div class=\"page\"><header class=\"hero\"><div class=\"brand\"><span class=\"brand-mark\" aria-hidden=\"true\">▦</span><span>rust-hic <span class=\"brand-sub\">/ run report</span></span></div><div class=\"hero-content\"><div><div class=\"kicker\">ANALYSIS REPORT <span>•</span> COMPLETE</div><h1>From reads<br>to contacts<span class=\"dot\">.</span></h1><p class=\"hero-copy\">A summary of the completed Hi-C processing run, from input through balanced contact matrices.</p></div><div class=\"hero-meta\"><div><span>ASSEMBLY</span><strong>");
    out.push_str(&escape(&r.assembly));
    out.push_str("</strong></div><div><span>INPUT</span><strong>");
    out.push_str(if r.input.mode == "fastq" {
        "Paired FASTQ"
    } else {
        "Pairs files"
    });
    out.push_str("</strong></div><div><span>ALIGNER</span><strong>");
    out.push_str(&escape(
        r.input.aligner.as_deref().unwrap_or("Not applicable"),
    ));
    out.push_str("</strong></div><div><span>ELAPSED</span><strong>");
    out.push_str(&duration(r.elapsed_seconds));
    out.push_str("</strong></div></div></div></header><main><section class=\"overview\" aria-label=\"Run overview\"><div class=\"section-heading\"><div><div class=\"kicker dark\">01 / OVERVIEW</div><h2>The run at a glance</h2></div><p>Counts are read pairs or contacts, not individual reads.</p></div><div class=\"cards\">");
    let (input_label, input_value, input_detail) = match r.counts.input_read_pairs {
        Some(n) => (
            "INPUT READ PAIRS",
            number(n),
            "Read pairs processed from FASTQ",
        ),
        None => (
            "INPUT PAIRS",
            number(r.counts.pairs_before_restriction),
            "Pairs entering downstream processing",
        ),
    };
    card(&mut out, input_label, &input_value, input_detail);
    card(
        &mut out,
        "FINAL CONTACTS",
        &number(r.counts.final_contacts),
        "After filtering and deduplication",
    );
    card(
        &mut out,
        "CIS CONTACTS",
        &percentage(r.counts.cis, r.counts.final_contacts),
        "Share of final contacts on the same chromosome",
    );
    card(
        &mut out,
        "DUPLICATES REMOVED",
        &percentage(r.counts.duplicates, r.counts.pairs_before_dedup),
        "Share of pairs entering pipeline deduplication",
    );
    out.push_str("</div>");
    let not_converged = r
        .resolutions
        .iter()
        .filter(|level| level.converged == Some(false))
        .count();
    let unusable = r
        .resolutions
        .iter()
        .filter(|level| level.converged.is_none())
        .count();
    if not_converged + unusable > 0 {
        let _ = write!(out, "<div class=\"callout alert\"><strong>Balance needs attention.</strong> {not_converged} resolution(s) did not converge; {unusable} have no usable bins. Review the matrix table before using balanced values.</div>");
    }
    out.push_str("</section><section class=\"panel\"><div class=\"section-heading\"><div><div class=\"kicker dark\">02 / CONTACTS</div><h2>Where the pairs went</h2></div><p>Counts from completed pipeline stages.</p></div><div class=\"table-wrap\"><table><thead><tr><th>Stage / category</th><th class=\"numeric\">Pairs</th><th class=\"numeric\">Share</th><th>Denominator or meaning</th></tr></thead><tbody>");
    if let Some(n) = r.counts.input_read_pairs {
        row(
            &mut out,
            "Input read pairs",
            n,
            None,
            "Paired reads processed from FASTQ",
        );
        row(
            &mut out,
            "Pairs entering processing",
            r.counts.pairs_before_restriction,
            Some(n),
            "Share of input read pairs; aligner selection may already have occurred",
        );
    } else {
        row(
            &mut out,
            "Pairs entering processing",
            r.counts.pairs_before_restriction,
            None,
            "Records entering restriction filtering or deduplication",
        );
    }
    if let Some(f) = &r.restriction.stats {
        row(
            &mut out,
            "Same-fragment pairs removed",
            f.same_fragment,
            Some(f.input_records),
            "Share of pairs entering restriction filtering",
        );
        row(
            &mut out,
            "After restriction filtering",
            f.retained_records,
            Some(f.input_records),
            "Share of pairs entering restriction filtering",
        );
    }
    row(
        &mut out,
        "Duplicates removed",
        r.counts.duplicates,
        Some(r.counts.pairs_before_dedup),
        "Share of pairs entering pipeline deduplication",
    );
    row(
        &mut out,
        "Final contacts",
        r.counts.final_contacts,
        Some(r.counts.pairs_before_dedup),
        "Share of pairs entering pipeline deduplication",
    );
    row(
        &mut out,
        "Cis",
        r.counts.cis,
        Some(r.counts.final_contacts),
        "Share of final contacts",
    );
    row(
        &mut out,
        "Cis ≥ 20 kb",
        r.counts.cis_20kb_plus,
        Some(r.counts.cis),
        "Share of cis contacts",
    );
    row(
        &mut out,
        "Cis < 20 kb",
        r.counts.cis_shorter_than_20kb,
        Some(r.counts.cis),
        "Share of cis contacts",
    );
    row(
        &mut out,
        "Trans",
        r.counts.trans,
        Some(r.counts.final_contacts),
        "Share of final contacts",
    );
    out.push_str("</tbody></table></div><p class=\"footnote\">Pairs entering processing are not a standard mapping rate: filtering may already have occurred within the aligner. Duplicate counts cover this pipeline's deduplication stage and may not include duplicates removed by an aligner.</p>");
    if let Some(f) = &r.restriction.stats {
        let _ = write!(out, "<div class=\"callout\"><strong>Restriction filtering enabled.</strong> Of {} same-fragment pairs removed, {} were classified as dangling, {} as self-ligation and {} as other orientations.</div>", number(f.same_fragment), number(f.dangling), number(f.self_ligation), number(f.unknown_mechanism));
    } else {
        out.push_str("<div class=\"callout neutral\"><strong>Restriction filtering was not enabled.</strong> No same-fragment pairs were removed by this stage.</div>");
    }
    out.push_str("</section>");
    if let Some(bins) = &r.counts.cis_distance_bins {
        out.push_str("<section class=\"panel\"><div class=\"section-heading\"><div><div class=\"kicker dark\">CIS DISTANCE</div><h2>Contacts by genomic distance</h2></div><p>Final cis contacts grouped by distance between ends.</p></div><div class=\"distance-chart\">");
        let max = bins.iter().map(|bin| bin.count).max().unwrap_or(0);
        for bin in bins {
            let label = bin.end_bp.map_or_else(
                || format!("≥ {}", resolution(bin.start_bp)),
                |end| format!("{} – {}", resolution(bin.start_bp), resolution(end)),
            );
            let width = if max == 0 {
                0.0
            } else {
                100.0 * bin.count as f64 / max as f64
            };
            let _ = write!(out, "<div class=\"distance-row\"><span>{}</span><div class=\"bar\" aria-hidden=\"true\"><span style=\"width:{width:.2}%\"></span></div><span class=\"numeric\">{}</span><span class=\"numeric subtle\">{}</span></div>", escape(&label), number(bin.count), percentage(bin.count, r.counts.cis));
        }
        out.push_str("</div><p class=\"footnote\">Intervals include the lower bound and exclude the upper bound; percentages use all final cis contacts. These raw counts are not normalized by genomic opportunity or bin width.</p></section>");
    } else {
        out.push_str("<section class=\"panel\"><div class=\"section-heading\"><div><div class=\"kicker dark\">CIS DISTANCE</div><h2>Contacts by genomic distance</h2></div></div><p class=\"footnote\">Distance, pair-type and strand distributions were not recorded by this run and cannot be reconstructed from its statistics files.</p></section>");
    }
    if let (Some(types), Some(strands)) = (&r.counts.pair_types, &r.counts.strand_orientations) {
        out.push_str("<section class=\"panel\"><div class=\"section-heading\"><div><div class=\"kicker dark\">CONTACT COMPOSITION</div><h2>Pair types and orientations</h2></div><p>Counts after deduplication.</p></div><div class=\"composition\"><div><h3>Pair type</h3><div class=\"table-wrap\"><table><thead><tr><th>Type</th><th class=\"numeric\">Contacts</th><th class=\"numeric\">Share</th></tr></thead><tbody>");
        for (kind, n) in types {
            let _ = write!(out, "<tr><th scope=\"row\">{}</th><td class=\"numeric\">{}</td><td class=\"numeric\">{}</td></tr>", escape(kind), number(*n), percentage(*n, r.counts.final_contacts));
        }
        out.push_str("</tbody></table></div></div><div><h3>Strand orientation</h3><div class=\"table-wrap\"><table><thead><tr><th>Strands</th><th class=\"numeric\">Contacts</th><th class=\"numeric\">Share</th></tr></thead><tbody>");
        for kind in ["++", "+-", "-+", "--", "other"] {
            if let Some(n) = strands.get(kind) {
                let _ = write!(out, "<tr><th scope=\"row\">{}</th><td class=\"numeric\">{}</td><td class=\"numeric\">{}</td></tr>", escape(kind), number(*n), percentage(*n, r.counts.final_contacts));
            }
        }
        out.push_str("</tbody></table></div></div></div><p class=\"footnote\">Strand and pair-type counts reflect the final pair order after flipping and sorting, for both cis and trans contacts.</p></section>");
    }
    out.push_str("<section class=\"panel\"><div class=\"section-heading\"><div><div class=\"kicker dark\">03 / MATRIX</div><h2>Resolution &amp; balance</h2></div><p>Check balance status before using normalized values.</p></div><div class=\"table-wrap\"><table><thead><tr><th>Resolution</th><th class=\"numeric\">Bins</th><th class=\"numeric\">Usable bins</th><th class=\"numeric\">Pixels</th><th>Balance status</th><th class=\"numeric\">Iterations</th></tr></thead><tbody>");
    for level in &r.resolutions {
        let (status, class) = balance_state(level);
        let iters = level.iterations.map_or("—".into(), |v| format!("{v:.0}"));
        let _ = write!(out, "<tr><th scope=\"row\">{}</th><td class=\"numeric\">{}</td><td class=\"numeric\">{} <span class=\"subtle\">({})</span></td><td class=\"numeric\">{}</td><td><span class=\"badge {}\">{}</span></td><td class=\"numeric\">{}</td></tr>", resolution(level.bp), number(level.bins as u64), number(level.usable_bins as u64), percentage(level.usable_bins as u64, level.bins as u64), number(level.pixels as u64), class, status, iters);
    }
    out.push_str("</tbody></table></div><p class=\"footnote\">Usable bins have a finite, positive balance weight. A resolution with no usable bins is reported separately even if its stored convergence flag is true.</p></section><section class=\"panel\"><div class=\"section-heading\"><div><div class=\"kicker dark\">04 / PERFORMANCE</div><h2>Time spent</h2></div><p>Wall-clock durations from the pipeline event log.</p></div><div class=\"timings\">");
    for t in &r.timings {
        let width = if r.elapsed_seconds > 0.0 {
            (t.seconds / r.elapsed_seconds * 100.0).clamp(0.0, 100.0)
        } else {
            0.0
        };
        let _ = write!(out, "<div class=\"time-row\"><div class=\"time-name\">{}</div><div class=\"bar\" aria-hidden=\"true\"><span style=\"width:{width:.2}%\"></span></div><div class=\"time-value\">{}</div></div>", escape(t.stage), duration(t.seconds));
    }
    let _ = write!(out, "</div><p class=\"footnote\">Total: {}. Alignment and chunk preparation overlap; their durations must not be added together. Total is measured up to report generation.</p></section>", duration(r.elapsed_seconds));
    out.push_str("<section class=\"panel\"><div class=\"section-heading\"><div><div class=\"kicker dark\">05 / RUN DETAILS</div><h2>Inputs and settings</h2></div><p>Paths and options used by this run.</p></div><dl class=\"details\">");
    let restriction = match r.restriction.mode {
        "arima_digest" => "Arima digest",
        "bed" => "Provided fragment BED",
        _ => "Disabled",
    };
    let entries = [
        ("Assembly", r.assembly.clone()),
        (
            "Aligner",
            r.input
                .aligner
                .clone()
                .unwrap_or_else(|| "Not applicable".into()),
        ),
        ("Restriction filtering", restriction.into()),
        (
            "Deduplication max mismatch",
            format!("{} bp", r.dedup_max_mismatch),
        ),
        ("Chromosome sizes", r.input.chromsizes.display().to_string()),
        (
            "Result directory",
            r.outputs.result_directory.display().to_string(),
        ),
    ];
    for (label, value) in entries {
        let _ = write!(out, "<dt>{}</dt><dd>{}</dd>", escape(label), escape(&value));
    }
    if let Some(reference) = &r.input.reference {
        let _ = write!(
            out,
            "<dt>Reference FASTA</dt><dd>{}</dd>",
            escape(&reference.display().to_string())
        );
    }
    if let Some(source) = &r.restriction.source {
        let _ = write!(
            out,
            "<dt>Restriction source</dt><dd>{}</dd>",
            escape(&source.display().to_string())
        );
    }
    if let Some(reused) = r.input.fastq_index_reused {
        let _ = write!(
            out,
            "<dt>FASTQ index</dt><dd>{}</dd>",
            if reused {
                "Reused"
            } else {
                "Built for this run"
            }
        );
    }
    for (i, path) in r.input.paths.iter().enumerate() {
        let _ = write!(
            out,
            "<dt>Input {}</dt><dd>{}</dd>",
            i + 1,
            escape(&path.display().to_string())
        );
    }
    out.push_str(
        "</dl><div class=\"output-line\"><strong>Output files</strong><span>contacts.mcool</span>",
    );
    if r.outputs.pairs.is_some() {
        out.push_str("<span>filtered-nodups.pairs</span>");
    }
    if r.outputs.bam.is_some() {
        out.push_str("<span>alignments.bam</span>");
    }
    if r.outputs.compressed_pairs.is_some() {
        out.push_str("<span>filtered-nodups.pairs.gz + .px2</span>");
    }
    out.push_str("<span>report.json</span></div>");
    if !r.chunks.is_empty() {
        let _ = write!(out, "<details class=\"chunks\"><summary>View {} alignment chunks</summary><div class=\"table-wrap\"><table><thead><tr><th>Chunk</th><th class=\"numeric\">Input read pairs</th><th class=\"numeric\">Selected pairs</th><th class=\"numeric\">Elapsed</th></tr></thead><tbody>", r.chunks.len());
        for c in &r.chunks {
            let _ = write!(out, "<tr><th scope=\"row\">{}</th><td class=\"numeric\">{}</td><td class=\"numeric\">{}</td><td class=\"numeric\">{}</td></tr>", c.chunk, number(c.input_pairs), c.selected_pairs.map_or("—".into(), number), duration(c.elapsed_seconds));
        }
        out.push_str("</tbody></table></div><p class=\"footnote\">Selected pairs are shown when reported by the aligner path. Chunk durations overlap.</p></details>");
    }
    out.push_str("</section><footer><span>rust-hic <span class=\"subtle\">/ analysis report</span></span><span>Generated by rust-hic ");
    out.push_str(escape(r.report_generator_version).as_str());
    out.push_str(" · Metrics describe this run; they do not by themselves establish data quality.</span></footer></main></div></body></html>\n");
    out
}

#[cfg(test)]
mod tests {
    use super::escape;

    #[test]
    fn escapes_user_controlled_html() {
        assert_eq!(escape("<&>\"'"), "&lt;&amp;&gt;&quot;&#39;");
    }
}
