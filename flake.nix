{
    description = "A very basic flake";

    inputs = {
        nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    };

    outputs = { self, nixpkgs }: 
    let 
        system = "x86_64-linux";
        starpuOverlay = f: p: {
            StarPU = p.callPackage ./starpu.nix {};
        };
        pkgs = import nixpkgs {
            inherit system;
            overlays = [ starpuOverlay ];
        };
    in
    {
        devShell = pkgs.mkShell {
            buildInputs = with pkgs; [
                StarPU
            ];
        };
    };
}
