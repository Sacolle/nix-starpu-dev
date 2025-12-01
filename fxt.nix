# from https://github.com/PhoqueEberlue/nixpkgs/blob/add-starpu/pkgs/by-name/st/starpu/package.nix
# TODO: adicionar opções de compilação extra, como debug
# https://files.inria.fr/starpu/doc/starpu.pdf
{ 
      # derivation dependencies
    lib,
    fetchurl,
    stdenv,
    autoreconfHook,


    perl,
    help2man,

    static ? false
}:
stdenv.mkDerivation (finalAttrs: {
    pname = "fxt";
    system = "x86_64-linux";
    version = "0.3.14";

    inherit static;

    src = fetchurl {
        url = "https://download.savannah.gnu.org/releases/fkt/fxt-${finalAttrs.version}.tar.gz";
        hash = "sha256-MX2Nkxdc2fJ+xDuDkLbSncZhFPBqp08jKYR9Sbqq6/I=";
    };

    nativeBuildInputs = [
        perl
        help2man
 #       autoreconfHook
    ];

    buildInputs = [
        perl
        help2man
  #      autoreconfHook
    ];

    configureFlags = (lib.optional finalAttrs.static ["CFLAGS=-fPIC" "--enable-static=yes" "--enable-shared=no"]); 


      # No need to add flags for CUDA, it should be detected by ./configure

      postConfigure = ''
        # Patch shebangs recursively because a lot of scripts are used
        shopt -s globstar
        patchShebangs --build **/*.sh
      '';

      postPatch = ''
        # Patch shebangs recursively because a lot of scripts are used
        shopt -s globstar
        patchShebangs --build **/*.sh
      '';

#      enableParallelBuilding = true;
      doCheck = true;
})



