SDK=/usr/lib/dart
#SDK=/stuff/install/eclipse-mars/dart-sdk/
PUB=$(SDK)/bin/pub 
DART=$(SDK)/bin/dart
DDC=$(SDK)/bin/dartdevc

release:
	${PUB} -v build --mode release

debug:
	${PUB} -v build --mode debug

upgrade:
	${PUB} upgrade

test:
	cd web; ./grid_test_suite.sh 1

test_long:
	cd web; ./grid_test_suite.sh 100

get:
	${PUB} get


webserver_ddc:
	echo http://localhost:8080/microcosm.html
	$(PUB) serve web/ --web-compiler=dartdevc
