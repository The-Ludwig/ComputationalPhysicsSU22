EXECUTABLE=alpha

all: build/exercise.zip

build/$(EXECUTABLE): $(EXECUTABLE).cpp build/Makefile | build 
	cmake --build build

build/report.pdf: build/output/pdf/barriers_test_100-Po-212-pot.pdf report.md | build
	pandoc report.md --number-sections --citeproc --bibliography=references.bib -f markdown -o build/report.pdf  

build/output/pdf/barriers_test_100-Po-212-pot.pdf: visualize.py build/output/barriers_test_100-Po-212-pot.npy | build
	python visualize.py

build/output/barriers_test_100-Po-212-pot.npy: build/$(EXECUTABLE) config.yaml | build
	build/$(EXECUTABLE)

build/exercise.zip: build/report.pdf | build
	zip -r build/exercise build/report.pdf 

build/conan.lock: conanfile.txt | build
	conan install -of build -if build .

build/Makefile: CMakeLists.txt build/conan.lock | build
	cmake -S . -B build
	
build: 
	mkdir -p build/output/pdf

# clean
clean: FORCE
	rm -rf build

FORCE:

.PHONY : all build FORCE clean cmake conan
