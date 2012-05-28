TOP=../..
PKG=${shell cd ${TOP};pwd}
SUITE=doRUnit.R
##PKGTEST=../../../PackageBuilds/SNPchip
##PKGCHECK=${shell cd ${PKGTEST};pwd}
#R=${R_HOME}/bin/R

all: inst test

inst: # Install package

	cd ${TOP}/..;\
	R CMD INSTALL --no-lock ${PKG}

test: # Run unit tests
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R --vanilla --slave < ${SUITE}

##check: # Run package checker
##	rm -rfv ${PKGTEST}
##	##mkdir ${PKGTEST}/${PKG}
##	rsync -Cavuzb --exclude '*.aux' --exclude '.DS*'  --exclude '*.dvi' --exclude '*.blg' --exclude '*.bbl' --exclude 'Rplots.pdf' --exclude '*-*.pdf' --exclude '*.tex' --exclude '*.so' --exclude 'scripts/test*' --exclude '.git*' --exclude 'push' --exclude '.auto*' ${PKG}/ ${PKGCHECK}
##	find . -name '*~' -exec rm {} \;
##	find . -name '#*' -exec rm {} \;
##	find . -name '.#*' -exec rm {} \;
##	find . -name '.Rhistory' -exec rm {} \;
##	###cd ${PKGTEST}
##	R CMD build --no-vignettes ${PKGCHECK}  ## fastest
##	R CMD check ${PKGCHECK}.tar.gz
