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
		    texstudio
		    python3Packages.pygments
                    (texlive.combine { inherit (texlive) scheme-medium 
		    disser 
		    collection-langcyrillic
		    collection-latexextra
		    collection-bibtexextra
		    babel-russian; })
                  ];


    QT_QPA_PLATFORM_PLUGIN_PATH="${pkgs.qt514.qtbase.bin}/lib/qt-${pkgs.qt514.qtbase.version}/plugins";
}
