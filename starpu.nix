# from https://github.com/PhoqueEberlue/nixpkgs/blob/add-starpu/pkgs/by-name/st/starpu/package.nix
# TODO: adicionar opções de compilação extra, como debug
# https://files.inria.fr/starpu/doc/starpu.pdf
{ 
  # derivation dependencies
  lib,
  fetchurl,
  stdenv,
  writableTmpDirAsHomeHook,
  autoreconfHook,
  # starpu dependencies
  hwloc,
  libuuid,
  libX11,
  fftw,
  fftwFloat, # Same than previous but with float precision
  pkg-config,
  libtool,
  simgrid,
  mpi,
  cudaPackages,

  python313,
  fxt,

  # Options
  buildMode ? "debug",
  maxBuffers ? 8,
  enableSimgrid ? false,
  enableMPI ? false,
  enableCUDA ? false,
  enableTrace ? false,
}:
stdenv.mkDerivation (finalAttrs: {
    pname = "StarPU";
    system = "x86_64-linux";
    version = "1.4.7";

    inherit buildMode;
    inherit maxBuffers;
    inherit enableSimgrid;
    inherit enableMPI;
    inherit enableCUDA;
    inherit enableTrace;

    src = fetchurl {
        url = "http://files.inria.fr/starpu/starpu-${finalAttrs.version}/starpu-${finalAttrs.version}.tar.gz";
        hash = "sha256-HrPfVRCJFT/m4LFyrZURhDS0qB6p6qWiw4cl0NtTsT4=";
    };
    nativeBuildInputs = [
        pkg-config
        hwloc
        libtool
        writableTmpDirAsHomeHook
        autoreconfHook

        python313
        fxt
    ] 
        ++ lib.optional finalAttrs.enableSimgrid simgrid
        ++ lib.optional finalAttrs.enableMPI mpi
        ++ lib.optional finalAttrs.enableCUDA cudaPackages.cudatoolkit;

    buildInputs = [
        libuuid
        libX11
        fftw
        fftwFloat
        hwloc

        python313
        fxt
    ]
        ++ lib.optional finalAttrs.enableSimgrid simgrid
        ++ lib.optional finalAttrs.enableMPI mpi
        ++ lib.optional finalAttrs.enableCUDA cudaPackages.cudatoolkit;

    

    configureFlags = [
        (lib.enableFeature true "quick-check")
        (lib.enableFeature false "build-examples")
        (lib.enableFeature finalAttrs.enableSimgrid "simgrid")

        (lib.enableFeature false "starpupy")

         # Static linking is mandatory for smpi
        (lib.enableFeature finalAttrs.enableMPI "mpi")
        (lib.enableFeature finalAttrs.enableMPI "mpi-check")
        (lib.enableFeature (!finalAttrs.enableMPI) "shared") 

        # (lib.optional finalAttrs.enableTrace "--prefix=${fxt}")
    ] ++ (
        if finalAttrs.buildMode == "debug" then 
            [ "--enable-debug" "--enable-verbose" "--enable-spinlock-check" ]
        else
            if finalAttrs.buildMode == "release" then 
                [ "--enable-fast" ]
            else builtins.throw "unrecognized build mode"
    ) ++ (lib.optional (finalAttrs.maxBuffers != 8) "--enable-maxbuffers=${toString finalAttrs.maxBuffers}");


      # No need to add flags for CUDA, it should be detected by ./configure

      patchPhase = ''
        # Patch shebangs recursively because a lot of scripts are used
        shopt -s globstar
        patchShebangs --build **/*.in
      '';

      postConfigure = ''
        # Patch shebangs recursively because a lot of scripts are used
        shopt -s globstar
        patchShebangs --build **/*.sh
      '';

      enableParallelBuilding = true;
      doCheck = true;
})



