all: _imagesc.so _wview.so plotn _speech.so _mojave.so

_imagesc.so: _imagesc.c
	gcc -O3 _imagesc.c -o _imagesc.so -fPIC -shared -I/usr/include/SDL2 -lSDL2 -lm

_wview.so: _wview.c
	gcc -O3 _wview.c -o _wview.so -I/usr/include/SDL2 -lSDL2 -lm -fPIC -shared -pthread

plotn: plotn.c
	gcc -O3 plotn.c -o plotn -I/usr/include/SDL -lSDL -lm

_speech.so: _speech.c
	gcc -O3 _speech.c -o _speech.so -fPIC -shared -lm

_mojave.so: _mojave.c
	gcc -O3 _mojave.c -o _mojave.so -fPIC -shared -I/usr/include/SDL2 -lSDL2 -lm -lSDL2_ttf


clean:
	rm _imagesc.so _wview.so plotn
