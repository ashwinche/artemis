
EIGEN=/usr/local/include/eigen3/
LBFGSpp=/home/ashwin/repo/libra/lib/LBFGSpp/include
LBFGS=$(BUILD_DIR)/../lib/lbfgs
LIBXML2=/usr/include/libxml2/
GLAD=/home/ashwin/repo/voxels/lib/glad/include

BUILD_DIR = $(shell pwd)
SRC_DIR = $(BUILD_DIR)/../src

GLFW = $(BUILD_DIR)/../lib/glfw-3.3.2

CC=gcc -O4
CXX=g++ -g -O4
RM=rm -f

IPATH    = -I$(VULKAN_SDK)/include -I$(GLFW)/include -I${CURDIR}/../include -I${EIGEN} -I$(LBFGSpp) -I$(LBFGS)/include -I$(LIBXML2) -I$(GLAD)
RPATH    = -Wl,-rpath,$(GLFW)/src
CPPFLAGS = -g $(IPATH)
LDFLAGS  = -g -L$(GLFW)/src 
LDLIBS   =
LIBS     = -lteem -lsfml-graphics -lsfml-window -lsfml-system -lglfw3 -lGL -lX11 -lXi -lXrandr -lXxf86vm -lXinerama -lXcursor -lrt -lm -pthread -ldl -lpng -lz -lxml2
CLIBS    = -lm

SRCS=$(wildcard $(SRC_DIR)/[!main]*.cpp)
DEPS=$(wildcard $(SRC_DIR)/*.d)
CSRCS=$(wildcard $(SRC_DIR)/*.c)
OBJS=$(subst .cpp,.o,$(SRCS))

# print ${EIGEN}
$(info    $(IPATH))

# dbg:
# 	ls $(EIGEN)
all: game cliplap viewnrrd tracking

game: $(OBJS) $(SRC_DIR)/main_game.o lbfgs
	$(CXX) $(CPPFLAGS) $(IPATH) $(LDFLAGS) -o $@ $(LBFGS)/lbfgs.o $(OBJS) $(SRC_DIR)/main_$@.o $(LDLIBS) $(RPATH) $(LIBS)

dconv: $(OBJS) $(SRC_DIR)/main_dconv.o lbfgs
	$(CXX) $(CPPFLAGS) $(IPATH) $(LDFLAGS) -o $@ $(LBFGS)/lbfgs.o $(OBJS) $(SRC_DIR)/main_$@.o $(LDLIBS) $(RPATH) $(LIBS)

cliplap: $(OBJS) $(SRC_DIR)/main_cliplap.o lbfgs
	$(CXX) $(CPPFLAGS) $(IPATH) $(LDFLAGS) -o $@ $(LBFGS)/lbfgs.o $(OBJS) $(SRC_DIR)/main_$@.o $(LDLIBS) $(RPATH) $(LIBS)

viewnrrd: $(OBJS) $(SRC_DIR)/main_viewnrrd.o lbfgs
	$(CXX) $(CPPFLAGS) $(IPATH) $(LDFLAGS) -o $@ $(LBFGS)/lbfgs.o $(OBJS) $(SRC_DIR)/main_$@.o $(LDLIBS) $(RPATH) $(LIBS)

tracking: $(OBJS) $(SRC_DIR)/main_tracking.o lbfgs
	$(CXX) $(CPPFLAGS) $(IPATH) $(LDFLAGS) -o $@ $(LBFGS)/lbfgs.o $(OBJS) $(SRC_DIR)/main_$@.o $(LDLIBS) $(RPATH) $(LIBS)

# game: $(OBJS) $(SRC_DIR)/main_$@.o lbfgs
# 	$(CXX) $(CPPFLAGS) $(IPATH) $(LDFLAGS) -o game $(LBFGS)/lbfgs.o $(OBJS) $(SRC_DIR)/main_$@.o $(LDLIBS) $(RPATH) $(LIBS)

lbfgs: $(LBFGS)/lbfgs.o

$(LBFGS)/lbfgs.o: $(LBFGS)/lbfgs.c
	# HELLO!
	$(CC) -c -I$(LBFGS)/include $< -o $(LBFGS)/lbfgs.o -lm 

%.o: %.cpp
	$(RM) $*.d
	$(CXX) $(IPATH) $(CPPFLAGS) -MM $*.cpp -MT $*.o > $*.d
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $*.cpp -o $*.o

# depend: .depend

# .depend: $(SRCS)
# 	# DEPEND
# 	$(RM) ./.depend
# 	$(CXX) $(IPATH) $(CPPFLAGS) -MM $^ -MT $(OBJS) >>./.depend;

clean:
	$(RM) $(SRC_DIR)/*.d $(SRC_DIR)/*.o game dconv cliplap viewnrrd tracking
	$(RM) $(OBJS)
	$(RM) $(LBFGS)/lbfgs.o

distclean: clean
	$(RM) *~ .depend

include $(DEPS)