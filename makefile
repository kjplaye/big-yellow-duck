all: _imagesc.so wview plotn _speech.so

_imagesc.so: _imagesc.c
	gcc _imagesc.c -o _imagesc.so -fPIC -shared -I/usr/include/SDL -lSDL -lm

wview: wview.c
	gcc wview.c -o wview -I/usr/include/SDL -lSDL -lm

plotn: plotn.c
	gcc plotn.c -o plotn -I/usr/include/SDL -lSDL -lm

_speech.so: _speech.c
	gcc _speech.c -o _speech.so -fPIC -shared -lm


clean:
	rm _imagesc.so wview plotn
