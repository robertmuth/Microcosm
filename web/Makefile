
SDK=/usr/lib/dart
DDC=$(SDK)/bin/dartdevc
PORT=8000
PUB=$(SDK)/bin/pub
DART=$(SDK)/bin/dart

release:
    ${PUB} build --mode release

debug:
    ${PUB} build --mode debug


test:
    cd web; ./grid_test_suite.sh 1

microcosm.html: debug.html
	./htmlpp.py <debug.html  >$@
