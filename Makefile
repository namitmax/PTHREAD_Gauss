CC=g++ -pthread

WARNINGS= #-mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

FLAGS=  #-g -pg -fsanitize=thread
FLAGS_OPTIMIZE= -O3 -ffast-math

all: a.out

a.out:  main.o input_output_tools.o calculations_tools.o
	$(CC) $(FLAGS) $(FLAGS_OPTIMIZE) $(WARNINGS) -o a.out main.o input_output_tools.o calculations_tools.o

input_output_tools.o: input_output_tools.cpp input_output_tools.h
	$(CC) $(FLAGS) $(FLAGS_OPTIMIZE) $(WARNINGS) -c input_output_tools.cpp

calculations_tools.o: calculations_tools.cpp calculations_tools.h input_output_tools.h
	$(CC) $(FLAGS) $(FLAGS_OPTIMIZE) $(WARNINGS) -c calculations_tools.cpp

main.o: main.cpp input_output_tools.h calculations_tools.h
	$(CC) $(FLAGS) $(FLAGS_OPTIMIZE) $(WARNINGS) -c  main.cpp

clean:
	rm -rf *.o prog 3 a.out
