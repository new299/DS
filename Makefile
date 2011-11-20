
ds_isolate: ds_isolate.cpp Image.h
	g++ ds_isolate.cpp -o ds_isolate -O3 -ltiff

ds_orientate: ds_orientate.cpp Image.h
	g++ ds_orientate.cpp -o ds_orientate -O3 -ltiff

ds_clean: ds_clean.cpp Image.h
	g++ ds_clean.cpp -o ds_clean -O3 -ltiff -fopenmp

ds_align: ds_align.cpp Image.h edit_distance.h
	#g++ ds_align.cpp -o ds_align -O3 -fopenmp
	icpc ds_align.cpp -o ds_align -O3 -openmp

ds_merge: ds_merge.cpp Image.h edit_distance.h
	g++ ds_merge.cpp -o ds_merge -O3 -ltiff -fopenmp

ds_merge_test: ds_merge_test.cpp Image.h edit_distance.h
	g++ ds_merge_test.cpp -o ds_merge_test -O3 -ltiff -fopenmp


ds_qc: ds_qc.cpp Image.h
	g++ ds_qc.cpp font.c spf.c -o ds_qc -g -ltiff -lSDLmain -lSDL -Wl,-framework,Cocoa -fopenmp
