CC = g++
CXXFLAGS += -std=c++11 -w
CXXFLAGS += -O3 -DNDEBUG -fopenmp
CXXOPENMP = 
LDFLAGS +=
LIBS += -lrt
# Specify the path of glib-core folder in SNAP source code files
GLIB = ../../Snap-6.0/glib-core
# Specify the path of snap-core folder in SNAP source code files
SNAP_CORE = ../../Snap-6.0/snap-core
# Specify the folder of phmap.h (the source code can be downloaded in https://github.com/greg7mdp/parallel-hashmap/)
PARALLEL_HASHMAP = ../../parallel-hashmap/parallel_hashmap

# Install SNAP library first, and specify the path of library if necessary
rKACS_Query: ./src/*.cpp
	$(CC) $(CXXFLAGS) -o rKACS_Query ./src/*.cpp -I$(SNAP_CORE) -I$(GLIB) -I$(PARALLEL_HASHMAP) -lsnap $(LDFLAGS) $(LIBS)
