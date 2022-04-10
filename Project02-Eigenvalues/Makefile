EXECUTABLE=$(shell pwd |  sed -e 's/.*\///' -e 's/\(.*\)/\L\1/' -e 's/project[0-9][0-9]-//')
PLOTS=disturbed_big.pdf
OUTPUT=disturbed_big.npy

# 'Debug' or 'Release'
BUILD_TYPE=Release

GCC_PATH=$(shell which g++)

PLOT_PATHS=$(addprefix build/plots/,$(PLOTS))
OUTPUT_PATHS=$(addprefix build/output/,$(OUTPUT))

$(info NAME=$(EXECUTABLE))
$(info OUTPUT_PATHS=$(OUTPUT_PATHS))
$(info PLOT_PATHS=$(PLOT_PATHS))
$(info g++=$(GCC_PATH))

.PHONY: all build FORCE clean

all: build/exercise.zip

build/exercise.zip: build/report.pdf | build
	cp build/report.pdf report.pdf
	zip -r build/exercise . --exclude build/\*
	rm report.pdf

build/report.pdf: report.md $(PLOT_PATHS) | build
	pandoc report.md --number-sections --citeproc --bibliography=references.bib --pdf-engine-opt=--shell-escape -f markdown -o build/report.pdf  

$(PLOT_PATHS)&: visualize.py $(OUTPUT_PATHS) | build
	python visualize.py

$(OUTPUT_PATHS)&: build/$(EXECUTABLE) config.yaml | build
	build/$(EXECUTABLE)

build/$(EXECUTABLE): build/Makefile $(EXECUTABLE).cpp | build 
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