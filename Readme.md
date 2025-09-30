# umwm-1d_rust

A Rust implementation of the 1D University of Miami Wave Model, translated from [Milan Curcic's original umwm-1d](https://github.com/umwm/umwm-1d). The model and source functions are largely based on the [University of Miami Wave Model](https://github.com/umwm/umwm) (UMWM), but simplified for one-dimensional applications.

## Prerequisites

This assumes familiarity with Rust. If you're new to Rust, [The Rust Programming Language](https://rust-book.cs.brown.edu/) by Brown University is an excellent starting resource.

## Usage

You can work with umwm-1d in two ways:

### Compiled Version

Run the main module directly from the root directory:
```bash
cargo run --release
```

Additional compilation options are available for different configurations.

### Notebook Version

You can use your environment to run rust in jupyter notebooks (that may be installed through your conda environment - say your_env)

To use Rust in Jupyter notebooks (don't need to be in a cond# umwm-1d_rust

A Rust implementation of the 1D University of Miami Wave Model, translated from [Milan Curcic's original umwm-1d](https://github.com/umwm/umwm-1d). The model and source functions are largely based on the [University of Miami Wave Model](https://github.com/umwm/umwm) (UMWM), but simplified for one-dimensional applications.

## Getting Started

Clone the repository:
```bash
git clone https://github.com/yourusername/umwm-1d_rust.git
cd umwm-1d_rust
```

## Prerequisites

This assumes familiarity with Rust. If you're new to Rust, [The Rust Programming Language](https://doc.rust-lang.org/book/) is an excellent starting resource.

## Usage

You can work with umwm-1d in two ways:

### Compiled Version

Run the main module directly from the root directory:
```bash
cargo run --no-default-features
```

Additional compilation options are available for different configurations.

### Notebook Version

You can run Rust in Jupyter notebooks using your existing environment (such as a conda environment).

First, install the Rust Jupyter kernel (works outside conda environments):
```bash
cargo install evcxr_jupyter
```

Then link the kernel to your conda environment and start Jupyter:
```bash
conda activate your_env
evcxr_jupyter --install
jupyter lab
```

**Note:** There are example Jupyter notebooks in the `notebooks` directory. Local library imports are not currently supported in the notebooks.

## Development Notes

While Rust can be verbose, the compilation process catches most errors upfront. After initial compilation, remaining bugs were primarily related to model formulation rather than memory management issues—a key advantage of Rust's ownership system.

## Disclaimer

This software is provided "as is" without warranty of any kind. The authors make no guarantees that it will work for your specific use case, environment, or expectations. Use at your own risk and discretion. The authors have no obligation to provide support, fix bugs, or maintain this code.a env)
```bash
cargo install evcxr_jupyter
```

To link the kernal with the notebook in your conda env
```
conda activate your_env
evcxr_jupyter --install
jupyter lab
```

**Note:** Local library imports are not currently supported in the notebooks directory.

## Development Notes

While Rust can be verbose, the compilation process catches most errors upfront. After initial compilation, remaining bugs were primarily related to model formulation rather than memory management issues—a key advantage of Rust's ownership system.
