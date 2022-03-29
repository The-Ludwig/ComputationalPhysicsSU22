EXECUTABLE=alpha

all: exercise.zip

build/$(EXECUTABLE): $(EXECUTABLE).cpp | conan cmake build 
	cmake --build build

report.pdf: | build
	pandoc report.md -o build/report.pdf

exercise.zip: report.pdf | build
	zip -r build/exercise build/report.pdf 

conan: conanfile.txt | build
	conan install -of build -if build .

cmake: CMakeLists.txt | build 
	cmake -S . -B build
	
build: 
	mkdir -p build

# clean
clean:
	rm -rf build

.PHONY : build clean cmake conan
