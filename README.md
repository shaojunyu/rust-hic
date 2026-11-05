<p align="center">
  <img src="docs/images/logo.png" alt="rust-hic logo" width="300">
</p>

<h1 align="center">rust-hic</h1>

`rust-hic` processes paired Hi-C FASTQ reads or pairs files into cleaned
contacts and a multi-resolution `.mcool` matrix. It includes BWA-MEM3 and also
supports external Chromap, BWA-MEM and BWA-MEM2. The pipeline is experimental.
Runs produce an offline HTML report and machine-readable statistics.

## Download and run

Download and extract the latest release:

```sh
curl -fL -o rust-hic-linux-x86_64.tar.gz https://github.com/shaojunyu/rust-hic/releases/latest/download/rust-hic-linux-x86_64.tar.gz
tar -xzf rust-hic-linux-x86_64.tar.gz
./rust-hic --help
./rust-hic example-config --minimal > config.yaml
# Edit config.yaml to point to your input files and output directory.
./rust-hic run config.yaml --check
./rust-hic run config.yaml
```

FASTQ input requires a reference FASTA, chromosome sizes and an existing
BWA-MEM2-compatible reference index. See [reference preparation](docs/usage.md#prepare-a-reference)
for details. To generate an annotated configuration with all YAML settings, run
`./rust-hic example-config > config.yaml`.

The `run` command defaults to local execution; `--mode slurm` submits the
workflow to a Slurm cluster. By default, the result includes `contacts.mcool`
and reports; BAM and pairs files are not retained.

The release executable targets x86-64 Linux; its built-in BWA-MEM3 requires
SSE4.2. External aligners must be installed separately. The archive also contains
license notices and documentation.

Restriction-fragment filtering is optional and disabled in the example by
default. See [restriction-fragment filtering](docs/usage.md#restriction-fragment-filtering)
before enabling it for your experiment.

- [Usage, configuration and commands](docs/usage.md)
- [Running across Slurm nodes](docs/slurm.md)
- [Building from source and development](docs/building.md)

## License

Project code is licensed under the [MIT License](LICENSE). Third-party
components retain their own licenses and copyright notices.
