{
  description = "Minimal Nix flake for building SFINCS version3";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs { inherit system; };
      mpi = pkgs.openmpi;
      mpiDev = pkgs.openmpi.dev;
      scalapack = pkgs.scalapack;
      scotch = pkgs.scotch.overrideAttrs (old: rec {
        version = "7.0.3";
        src = pkgs.fetchFromGitLab {
          domain = "gitlab.inria.fr";
          owner = "scotch";
          repo = "scotch";
          rev = "v${version}";
          sha256 = "0krs485gsf6vlhgkh2w96mspbmysmsalj4cnnjp9xk63w39qgcy2";
        };
      });
      mumps = (pkgs.mumps.override {
        inherit scotch;
        mpiSupport = true;
        withPtScotch = true;
      }).overrideAttrs (old: rec {
        version = "5.6.2";
        src = pkgs.fetchzip {
          url = "https://mumps-solver.org/MUMPS_${version}.tar.gz";
          sha256 = "1dwx4virk4a53pk5ny3v88c3xa8d3ym7157x63cqyq5qky6hy67m";
        };
      });
      petsc = (pkgs.petsc.override {
        hdf5-support = true;
        petsc-optimized = true;
      }).overrideAttrs (old: rec {
        version = "3.21.1";
        src = pkgs.fetchzip {
          url = "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${version}.tar.gz";
          sha256 = "12nczdxk61vy7gr73p8xbfpwic504v713285xhhrb8qabw6135pr";
        };
        buildInputs = old.buildInputs ++ [
          mumps
          scalapack
          scotch
        ];
        configureFlags = old.configureFlags ++ [
          "--with-mumps-include=${mumps}/include"
          "--with-mumps-lib=[${mumps}/lib/libdmumps.so,${mumps}/lib/libmumps_common.so,${mumps}/lib/libpord.so]"
          "--with-scalapack-lib=${scalapack}/lib/libscalapack.so"
        ];
      });
      hdf5 = pkgs.hdf5-fortran;
      netcdf = pkgs.netcdf;
      netcdffortran = pkgs.netcdffortran;
      python = pkgs.python3;
      hdf5Tools = pkgs.hdf5;
      buildTools = [
        pkgs.gnumake
        pkgs.pkg-config
        mpiDev
        netcdf
        netcdffortran
      ];
      buildLibs = [
        mpi
        petsc
        mumps
        scalapack
        scotch
        hdf5
        hdf5.dev
        netcdf
        netcdffortran
        pkgs.curl
        pkgs.zlib
      ];
      commonEnv = ''
        export PETSC_DIR=${petsc}
        export PETSC_ARCH=
        export PATH=${mpiDev}/bin:${mpi}/bin:${netcdffortran}/bin:${netcdf}/bin:$PATH
        export PKG_CONFIG_PATH=${petsc}/lib/pkgconfig:${hdf5.dev}/lib/pkgconfig:$PKG_CONFIG_PATH
      '';
    in {
      devShells.${system}.default = pkgs.mkShell {
        packages = buildTools ++ buildLibs;
        shellHook = commonEnv;
      };

      packages.${system}.default = pkgs.stdenv.mkDerivation {
        pname = "sfincs-version3";
        version = "unstable";
        src = ./.;
        strictDeps = true;
        doCheck = true;

        nativeBuildInputs = buildTools;
        buildInputs = buildLibs;
        nativeCheckInputs = [
          python
          hdf5Tools
        ];
        configurePhase = ":";

        buildPhase = ''
          ${commonEnv}
          export SFINCS_SYSTEM=nix
          make -C fortran/version3 clean
          make -C fortran/version3 -j$NIX_BUILD_CORES
        '';

        checkPhase = ''
          ${commonEnv}
          export SFINCS_SYSTEM=nix
          patchShebangs fortran/version3/examples
          find fortran/version3/examples -name 'job.CI' -execdir ln -sf job.CI job.nix ';'
          make -C fortran/version3 test
        '';

        installPhase = ''
          mkdir -p $out/bin
          cp fortran/version3/sfincs $out/bin/
        '';
      };
    };
}
