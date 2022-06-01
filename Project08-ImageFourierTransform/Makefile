EXECUTABLE=$(shell pwd |  sed -e 's/.*\///' -e 's/\(.*\)/\L\1/' -e 's/project[0-9][0-9]-//')
PLOTS=times_lin.pdf times_p2.pdf 
OUTPUT=times_lin.npy times_p2.npy

# 'Debug' or 'Release'
BUILD_TYPE=Release

GCC_PATH=$(shell which g++)

PLOT_PATHS=$(addprefix build/plots/,$(PLOTS))
OUTPUT_PATHS=$(addprefix build/output/,$(OUTPUT))

$(info NAME=$(EXECUTAEigen3::EigenBLE))
$(info OUTPUT_PATHS=$(OUTPUT_PATHS))
$(info PLOT_PATHS=$(PLOT_PATHS))
$(info g++=$(GCC_PATH))

.PHONY: all build FORCE clean

all: build/exercise.zip build/presentation.pdf

build/exercise.zip: build/report.pdf | build
	cp build/report.pdf report.pdf
	zip -r build/exercise . --exclude build/\*
	rm report.pdf


TEXOPTIONS = -lualatex \
		--output-directory=build \
		--interaction=nonstopmode \
		--halt-on-error \
		--use-make \
		--synctex=1 \

TEXPREFIX = TEXINPUTS=build: \
	BIBINPUTS=build: \
	max_print_line=1048576 \

# Use this to work on tex-document, it will be updated continuesly by latexmk
work: FORCE | build
	$(TEXPREFIX) latexmk -pvc $(TEXOPTIONS) report.tex

presentation: FORCE | build
	$(TEXPREFIX) latexmk -pvc $(TEXOPTIONS) presentation.tex

build/presentation.pdf: presentation.tex $(PLOT_PATHS) build/plots/hole_dft.png | build
	$(TEXPREFIX) latexmk $(TEXOPTIONS) presentation.tex

build/plots/hole_dft.png: build/img_dft | build
	./fft_all_images.bash

build/report.pdf: report.tex $(PLOT_PATHS) build/plots/hole_dft.png| build
	$(TEXPREFIX) latexmk $(TEXOPTIONS) report.tex

$(PLOT_PATHS)&: visualize.py $(OUTPUT_PATHS) | build
	python visualize.py

$(OUTPUT_PATHS)&: build/$(EXECUTABLE) | build
	build/$(EXECUTABLE)

build/$(EXECUTABLE) build/img_dft: build/Makefile FORCE | build
#	cmake -DNAME:STRING=$(EXECUTABLE) --build build 
	$(MAKE) -C build

build/Makefile: build/conan_paths.cmake CMakeLists.txt | build
	CXX=$(GCC_PATH) cmake -S . -B build -DNAME:STRING=$(EXECUTABLE) -DCMAKE_BUILD_TYPE:STRING=$(BUILD_TYPE)

build/conan%: conanfile.txt | build
	conan install -of build -if build .
	
build: FORCE
	mkdir -p build/output build/plots

# clean
clean: FORCE
	rm -rf build

FORCE: