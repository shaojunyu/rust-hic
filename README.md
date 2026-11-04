<p align="center">
  <img src="docs/images/logo.png" alt="rust-hic logo" width="240">
</p>

<h1 align="center">rust-hic</h1>

`rust-hic` processes paired Hi-C FASTQ reads or pairs files into cleaned
contacts and a multi-resolution `.mcool` matrix. It includes BWA-MEM3 and also
supports external Chromap, BWA-MEM and BWA-MEM2. The pipeline is experimental.

## Download and run

Download `rust-hic-linux-x86_64` from the latest [GitHub Release](../../releases/latest),
make it executable, and run it:

```sh
chmod +x rust-hic-linux-x86_64
mv rust-hic-linux-x86_64 rust-hic
./rust-hic --help
./rust-hic example-config > config.yaml
# Edit config.yaml to point to your input files and output directory.
./rust-hic run config.yaml --check
./rust-hic run config.yaml
```

The release executable targets x86-64 Linux; its built-in BWA-MEM3 requires
SSE4.2. External aligners must be installed separately. Download the release
archive for bundled license notices and documentation.

- [Usage, configuration and commands](docs/usage.md)
- [Building from source and development](docs/building.md)

## License

Project code is licensed under the [MIT License](LICENSE). Third-party
components retain their own licenses and copyright notices.
