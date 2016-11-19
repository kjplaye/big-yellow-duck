all: _imagesc.so _wview.so plotn _speech.so

_imagesc.so: _imagesc.c
	gcc _imagesc.c -o _imagesc.so -fPIC -shared -I/usr/include/SDL -lSDL -lm

_wview.so: _wview.c
	gcc _wview.c -o _wview.so -I/usr/include/SDL -lSDL -lm -fPIC -shared

plotn: plotn.c
	gcc plotn.c -o plotn -I/usr/include/SDL -lSDL -lm

_speech.so: _speech.c
	gcc _speech.c -o _speech.so -fPIC -shared -lm


clean:
	rm _imagesc.so _wview.so plotn
