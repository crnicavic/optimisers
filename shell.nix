with import <nixpkgs> {};
mkShell {
	buildInputs = [
		pkg-config
		gcc
		clang
		gdb
		cmake
		gnumake
		xorg.libXext
		xorg.libXrandr
		xorg.libX11
		xorg.libXcursor
		xorg.libXfixes
		libxkbcommon
		libdrm
		libGL
		libGLU
		SDL2
		SDL
		python3
		python312Packages.pip
		];
}
