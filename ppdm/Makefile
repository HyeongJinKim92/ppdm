CC = g++
BIN=PPDM
LIBRARIES=pthread
SOURCES_UTIL=util/*.cpp
OBJECTS_UTIL=util/*.o
SOURCES_NTL=NTL/*.cpp
OBJECTS_NTL=NTL/*.o
SOURCES_CIRCUIT=circuit/*.cpp
OBJECTS_CIRCUIT=circuit/*.o
SOURCES_PARALLEL=parallel_srange.cpp parallel_stopk.cpp parallel_sknn.cpp parallel_classification.cpp parallel_clustering.cpp parallel_protocol.cpp
INCLUDE=..

all: ${BIN}

PPDM: ${OBJECTS_NTL} ${OBJECTS_UTIL} ${OBJECTS_CIRCUIT} 
	${CC} -o PPDM ${CFLAGS} ${OBJECTS_NTL} ${OBJECTS_UTIL} ${OBJECTS_CIRCUIT} \
	-l${LIBRARIES} ${SOURCES_PARALLEL} ppdm.cpp paillier.c protocol.c sknn.c stopk.c srange.c G_sRange.c setting.c classification.c clustering.c main.c \
	-L/usr/include -lgmp -lm -Wno-write-strings

${OBJECTS_UTIL}: ${SOURCES_UTIL}$  
	@cd util; ${CC} -c -I${INCLUDE} ${CFLAGS} *.cpp

${OBJECTS_NTL}: ${SOURCES_NTL}$  
	@cd NTL; ${CC} -c -I${INCLUDE} ${CFLAGS} *.cpp

${OBJECTS_CIRCUIT}: ${SOURCES_CIRCUIT}$  
	@cd circuit; ${CC} -c -I${INCLUDE} ${CFLAGS} *.cpp

clean:
	rm -rf ${OBJECTS_UTIL} ${OBJECTS_NTL} ${OBJECTS_CIRCUIT} ${BIN}
