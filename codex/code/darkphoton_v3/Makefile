CXX = g++
CXXFLAGS = -O2 -std=c++17 $(shell root-config --cflags)
LDLIBS = $(shell root-config --libs)

.PHONY: all run reports
all: v3 report
v3: v3.cpp redo.cpp configuration.hxx physics_base.hxx algebra.hxx
	$(CXX) $(CXXFLAGS) v3.cpp $(LDLIBS) -o v3
report: report.cpp
	$(CXX) $(CXXFLAGS) report.cpp $(LDLIBS) -o report
output/danilov_page3.svg: ../../papers/darkphoton_v3/danilov2019.pdf
	pdftocairo -f 3 -l 3 -svg $< $@
run: all output/danilov_page3.svg
	./v3 config.conf > output/run_receipt.txt 2>&1
	./report
reports: run
	cd ../.. && TECTONIC_CACHE_DIR=code/darkphoton_v3/tex-cache XDG_CACHE_HOME=code/darkphoton_v3/tex-cache code/park_production_v2/tools/tectonic --only-cached --keep-logs report_darkphoton_v3_en.tex > code/darkphoton_v3/output/latex_en_receipt.txt 2>&1
	cd ../.. && TECTONIC_CACHE_DIR=code/darkphoton_v3/tex-cache XDG_CACHE_HOME=code/darkphoton_v3/tex-cache code/park_production_v2/tools/tectonic --only-cached --keep-logs report_darkphoton_v3_es.tex > code/darkphoton_v3/output/latex_es_receipt.txt 2>&1
	cp ../../report_darkphoton_v3_en.pdf output/
	cp ../../report_darkphoton_v3_es.pdf output/
