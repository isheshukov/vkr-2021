# Graduation thesis 2021

## Build
For dependency management you need to use [Nix](https://nixos.org/download.html).
Install Nix and inside the root of this repository run

```
nix-shell
mkdir build && cd build
cmake ..
make

```

To build release for Linux, Windows, and PDFs for graduation thesis text and slides run inside the root of this repository

```
make
```
Produced files will be in newly created `dist/` folder.

## Run

```
./build/app/tropical
```
