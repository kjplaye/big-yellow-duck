all: _imagesc.so _wview.so plotn _speech.so

_imagesc.so: _imagesc.c
	gcc -O3 _imagesc.c -o _imagesc.so -fPIC -shared -I/usr/include/SDL -lSDL -lm

_wview.so: _wview.c
	gcc -O3 _wview.c -o _wview.so -I/usr/include/SDL -lSDL -lm -fPIC -shared -pthread

plotn: plotn.c
	gcc -O3 plotn.c -o plotn -I/usr/include/SDL -lSDL -lm

_speech.so: _speech.c
	gcc -O3 _speech.c -o _speech.so -fPIC -shared -lm


clean:
	rm _imagesc.so _wview.so plotn
