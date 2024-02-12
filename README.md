# dsa: Deep Sequencing Analysis

## Getting Started

dsa is a command line application developed and tested on Microsoft Windows
11 and WSL2. Some of the code includes AVX2 SIMD instructions it requires
an x86 CPU (i.e., Intel or AMD CPUs but not Apple silicon).

## Prerequisites

The Windows version can be built with recent versions (>=17) of Visual Studio
Community Edition. Compilation on Linux requires a C++ compiler and standard 
library with support for C++20 features. The Linux build was tested in WSL2
with clang++ version 12. dsa has no external dependencies.

## Installing

For a Windows install, we recommend simply downloading the binary installer
(dsa-setup.msi). For the Linux build, clone the Git repository,
cd into it, run make, and copy the executable bin/dsa to a
directory in your shell's path. For example...

```
$cd ~
$git clone https://github.com/baileych-bi/dsa-win64
$cd dsa-win64
$make
$sudo cp bin/dsa /usr/local/bin
$dsa --help
```

## Usage

The dsa documentation available from running dsa with the --help option
is quite thorough. Examples of scripts for batch analysis using Powershell
and Bash are available in the 'example' folder along with sample .fastq files
(unzip after downloading) and sample output.

## Author

Charles C Bailey

## Acknowledgements

dsa uses a Windows implementation of getopt from the
mingw-w64 package for command line argument processing
on the Windows platform. See local-getopt.h for license
and authorship of that header library.
