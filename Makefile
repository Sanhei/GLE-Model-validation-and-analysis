# Makefile for program: Filament Dynamcis
# Benjamin Dalton: 18-09-2016
# need to  export PYTHONPATH="/home/physikqy/.conda/envs/clib/lib/python3.9/site-packages"
 
export PYTHONPATH="/home/physikqy/env/CM/lib/python3.9/site-packages"

EX_NAME = memoryextraction
CPP_FLAGS = -c -std=c++14 -O1

OBJS = main.o histogram_constructor.o correlation.o memory.o passage_time.o figure_plot.o

$(EX_NAME) : $(OBJS)
	@ export PYTHONPATH="/home/physikqy/.conda/envs/clib/lib/python3.9/site-packages"; \
	g++ -o $(EX_NAME) $(OBJS) -I/usr/include/python3.9 -I /home/physikqy/.conda/envs/clib/lib/python3.9/site-packages/numpy/core/include -L /home/physikqy/.conda/envs/clib/lib -lpython3.9 -lpthread -lutil -ldl -Xlinker -export-dynamic -lpython3.9 -lfftw3 -lfftw3_threads -lm -lstdc++fs  


main.o : main.cpp
	@ export PYTHONPATH="/home/physikqy/.conda/envs/clib/lib/python3.9/site-packages"; \
	g++ $(CPP_FLAGS) main.cpp -lstdc++fs


histogram_constructor.o : memory_extraction/histogram_constructor.cpp memory_extraction/potential.h
	g++ $(CPP_FLAGS) memory_extraction/histogram_constructor.cpp 

correlation.o : memory_extraction/correlation.cpp memory_extraction/correlation.h
	g++ $(CPP_FLAGS) memory_extraction/correlation.cpp -lfftw3 -lfftw3_threads -lm   
	
memory.o: memory_extraction/memory.cpp
	g++ $(CPP_FLAGS) memory_extraction/memory.cpp
	
passage_time.o : Passage_time/passage_time.cpp Passage_time/passage_time.h
	g++ $(CPP_FLAGS) Passage_time/passage_time.cpp	



figure_plot.o: figure_plot.cpp figure_plot.h 
	@ export PYTHONPATH="/home/physikqy/.conda/envs/clib/lib/python3.9/site-packages"; \
	g++ $(CPP_FLAGS) figure_plot.cpp -I/usr/include/python3.9 -I /home/physikqy/.conda/envs/clib/lib/python3.9/site-packages/numpy/core/include -L /home/physikqy/.conda/envs/clib/lib -lpython3.9 -lpthread -lutil -ldl -Xlinker -export-dynamic -lpython3.9 
clean:
	rm -f core $(EX_NAME) $(OBJS)
	
