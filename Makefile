CXX=clang++
CXXFLAGS=-g -std=c++17 -O3
LDFLAGS=-ligraph -lxml2

ifneq ($(CXX),icc)
LDFLAGS+=-lblas -lm
else
LDFLAGS+=${MKL_SHLIB}
endif

RM=rm -rf
MKDIR=mkdir -p

SRCS=$(wildcard src/*.cpp)
OBJS=$(subst src/,obj/,$(subst .cpp,.o,$(SRCS)))

OUT=simulation

.PHONY : all
all : bin/$(OUT)

.PHONY: clean
clean:
	@$(RM) obj/

.PHONY: distclean
distclean:
	@$(RM) obj/ bin/


bin/$(OUT): $(OBJS)
	@$(MKDIR) bin/
	$(CXX) -o bin/$(OUT) $^ $(LDFLAGS)

obj/%.o: src/%.cpp
	@$(MKDIR) obj/
	$(CXX) $(CXXFLAGS) -o $@ -c $<
