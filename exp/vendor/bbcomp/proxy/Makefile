
MODE = -O3 -DNDEBUG
#MODE = -g -DDEBUG

FLAGS = -I../library/ ${MODE} -Wall  -fPIC


proxy: proxy.cpp socket.cpp libbbcomp.so
	${CXX} ${FLAGS} proxy.cpp socket.cpp libbbcomp.so -ldl -o proxy

libbbcomp.so:
	cp ../library/libbbcomp.so .

clean:
	rm -f proxy libbbcomp.so
