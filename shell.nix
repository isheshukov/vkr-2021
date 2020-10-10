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
    ];

    buildInputs = with pkgs; [
                    mesa
		    mesa_glu
                    qt514.full
                    qtcreator
                  ];

    QT_QPA_PLATFORM_PLUGIN_PATH="${pkgs.qt514.qtbase.bin}/lib/qt-${pkgs.qt514.qtbase.version}/plugins";
}
