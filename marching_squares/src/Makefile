CC=g++
FLAGS=$$(pkg-config opencv4 --libs --cflags) -lstdc++ -lgdal

%: %.cpp
	$(CC) -o $@ $^ $(FLAGS) -pg
