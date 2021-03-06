let pkgs = import <nixpkgs> {};
in
pkgs.clangStdenv.mkDerivation {
    name = "vkr-shell";

    nativeBuildInputs = with pkgs; [
      gdb
      pkg-config
      cmake
      clang-analyzer
      clang-tools
      ginac
      cln
      valgrind
      zip
    ];

    buildInputs = with pkgs; [
                    ginac
		    cln
                    mesa
                    mesa_glu
                    qt514.full
                    qtcreator
		    dot2tex
		    texstudio
		    python3Packages.pygments
		    texlive.combined.scheme-full
                  ];


    QT_QPA_PLATFORM_PLUGIN_PATH="${pkgs.qt514.qtbase.bin}/lib/qt-${pkgs.qt514.qtbase.version}/plugins";
}
