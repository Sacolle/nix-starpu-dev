{
    description = "A very basic flake";

    inputs = {
        nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    };

    outputs = { self, nixpkgs }: 
    let 
        system = "x86_64-linux";
        starpuOverlay = f: p: {
            StarPU = p.callPackage ./starpu.nix { maxBuffers = 10; };
        };
        pkgs = import nixpkgs {
            inherit system;
            overlays = [ starpuOverlay ];
        };
    in
    {
        devShells.${system}.default = pkgs.mkShell {
            buildInputs = with pkgs; [
                pkg-config
                StarPU
                hwloc
                # for bash to work properlly inside vscode
                bashInteractive
                gdb
                valgrind
            ];
            shellHook = ''
              export SHELL=/run/current-system/sw/bin/bash
            '';
        };
    };
}
