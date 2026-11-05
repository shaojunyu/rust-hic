//! Optional preservation of aligner-emitted BAM records.
use anyhow::{bail, Context, Result};
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write as _};
use std::{
    fs::{File, OpenOptions},
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
};

pub(crate) type Writer = bam::io::Writer<noodles_bgzf::io::Writer<BufWriter<File>>>;

pub(crate) fn index_header(index: &bwa_mem3_rs::BwaIndex) -> Result<sam::Header> {
    let mut text = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    for (name, len) in index.contigs() {
        if name.is_empty() || name.contains(['\t', '\r', '\n']) || len <= 0 {
            bail!("invalid reference sequence in BWA index")
        }
        use std::fmt::Write as _;
        writeln!(text, "@SQ\tSN:{name}\tLN:{len}")?;
    }
    text.parse().context("construct BAM reference header")
}

fn create(path: &Path, header: &sam::Header) -> Result<Writer> {
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .with_context(|| format!("create BAM {}", path.display()))?;
    let mut writer = bam::io::Writer::new(BufWriter::with_capacity(1024 * 1024, file));
    writer.write_header(header).context("write BAM header")?;
    Ok(writer)
}

pub(crate) fn create_index_writer(path: &Path, index: &bwa_mem3_rs::BwaIndex) -> Result<Writer> {
    create(path, &index_header(index)?)
}

pub(crate) fn finish(writer: &mut Writer) -> Result<()> {
    writer.try_finish().context("finish BAM stream")?;
    // try_finish writes the BGZF EOF marker after flushing; flush the
    // underlying buffered file as well so late write errors are reported.
    writer.get_mut().flush().context("flush BAM file")?;
    Ok(())
}

/// Convert all SAM records emitted by an external aligner without applying pair filters.
pub(crate) fn sam_to_bam<R: BufRead>(input: R, output: &Path) -> Result<()> {
    let mut reader = sam::io::Reader::new(input);
    let header = reader.read_header().context("read aligner SAM header")?;
    if header.reference_sequences().is_empty() {
        bail!("aligner SAM header has no reference sequences")
    }
    let mut writer = create(output, &header)?;
    for record in reader.records() {
        let record = record.context("read aligner SAM record")?;
        writer
            .write_alignment_record(&header, &record)
            .with_context(|| format!("write aligner BAM record {:?}", record.name()))?;
    }
    finish(&mut writer).context("finish aligner BAM")?;
    Ok(())
}

/// Merge unsorted BAM shards, preserving their raw records and one reference dictionary.
pub(crate) fn merge(parts: &[PathBuf], output: &Path) -> Result<()> {
    if parts.is_empty() {
        bail!("no BAM chunks to merge")
    }
    let mut writer = None;
    let mut first_header: Option<sam::Header> = None;
    for part in parts {
        let file =
            File::open(part).with_context(|| format!("open BAM chunk {}", part.display()))?;
        let mut reader = bam::io::Reader::new(BufReader::with_capacity(1024 * 1024, file));
        let mut header = reader
            .read_header()
            .with_context(|| format!("read BAM header {}", part.display()))?;
        if let Some(expected) = &first_header {
            if !expected
                .reference_sequences()
                .iter()
                .eq(header.reference_sequences().iter())
            {
                bail!("BAM reference dictionaries differ in {}", part.display())
            }
        } else {
            use sam::header::record::value::map::header::tag;
            let hd = header.header_mut().get_or_insert_with(Default::default);
            hd.other_fields_mut()
                .insert(tag::SORT_ORDER, "unsorted".into());
            hd.other_fields_mut().shift_remove(&tag::SUBSORT_ORDER);
            writer = Some(create(output, &header)?);
            first_header = Some(header);
        }
        io::copy(
            reader.get_mut(),
            writer.as_mut().context("missing BAM writer")?.get_mut(),
        )
        .with_context(|| format!("copy BAM records from {}", part.display()))?;
    }
    finish(&mut writer.context("missing BAM writer")?).context("finish merged BAM")?;
    Ok(())
}

/// Duplicate a read stream into another consumer, propagating all I/O failures.
pub(crate) struct Tee<R, W> {
    pub input: R,
    pub copy: W,
}
impl<R: Read, W: Write> Read for Tee<R, W> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.input.read(buf)?;
        self.copy.write_all(&buf[..n])?;
        Ok(n)
    }
}
