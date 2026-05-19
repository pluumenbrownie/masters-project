{
  description = "Flake containing all needed dependencies to reproduce this masters thesis.";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    rust-overlay.url = "github:oxalica/rust-overlay";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = {
    nixpkgs,
    rust-overlay,
    flake-utils,
    ...
  }:
    flake-utils.lib.eachDefaultSystem (
      system: let
        overlays = [(import rust-overlay)];
        pkgs = import nixpkgs {
          inherit system overlays;
        };

        python = pkgs.python3.override {
          self = python;
          packageOverrides = pyfinal: pyprev: {
            # garage_admin_sdk = pyfinal.callPackage ./garage_sdk.nix {};
          };
        };
      in {
        devShells.default = with pkgs;
          mkShell {
            buildInputs = [
              (rust-bin.stable.latest.default.override {
                extensions = ["rust-src"];
              })
              bacon
              gcc
              cargo-flamegraph
              gnuplot

              (python.withPackages
                (pypkgs:
                  with pypkgs; [
                    torch
                    numpy
                  ]))
            ];

            shellHook = ''
              alias ls=eza
              alias find=fd
            '';
          };
      }
    );
}
