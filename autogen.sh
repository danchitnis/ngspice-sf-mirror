#!/bin/sh
# Configuration script for ngspice.
#
# This script performs initial configuration of ngspice source
# package.
#
#
# temp-adms.ac: modified configure.ac if --adms is selected
# for temporary use by autoconf, will be deleted automatically
# configure.ac stays untouched

PROJECT=ngspice

# Exit variable
DIE=0


# Check for Mac OS X
uname -a | grep -q "Darwin"
if [ $? -eq 0 ]; then
    LIBTOOLIZE=glibtoolize
else
    LIBTOOLIZE=libtoolize
fi

help()
{
    echo
    echo "$PROJECT autogen.sh help"
    echo
    echo "--help     -h: print this file"
    echo "--version  -v: print version"
    echo
}

version()
{
    echo
    echo "$PROJECT autogen.sh 1.0"
    echo
}

error_and_exit()
{
    echo "Error: $1"
    exit 1
}

check_awk()
{
    (awk --version) < /dev/null > /dev/null 2>&1 ||
    (awk -W version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have awk installed to compile $PROJECT with --adms."
	exit 1
    }
}

check_autoconf()
{
    (autoconf --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have autoconf installed to compile $PROJECT."
	echo "See http://www.gnu.org/software/automake/"
	echo "(newest stable release is recommended)"
	DIE=1
    }

    ($LIBTOOLIZE --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have libtool installed to compile $PROJECT."
	echo "See http://www.gnu.org/software/libtool/"
	echo "(newest stable release is recommended)"
	DIE=1
    }

    (automake --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have automake installed to compile $PROJECT."
	echo "See http://www.gnu.org/software/automake/"
	echo "(newest stable release is recommended)"
	DIE=1
    }
}


case "$1" in
    "--adms" | "-a")
        echo "Warning: adms is no longer available, ignored!"
        ;;

    "--help" | "-h")
        help
        exit 0
        ;;

    "--version" | "-v")
        version
        exit 0
        ;;

    *)
        ;;
esac


check_autoconf

if [ "$DIE" -eq 1 ]; then
    exit 1
fi

[ -f "DEVICES" ] || {
    echo "You must run this script in the top-level $PROJECT directory"
    exit 1
}


echo "Running $LIBTOOLIZE"
$LIBTOOLIZE --copy --force \
    || error_and_exit "$LIBTOOLIZE failed"

echo "Running aclocal $ACLOCAL_FLAGS"
aclocal $ACLOCAL_FLAGS --force -I m4 \
    || error_and_exit "aclocal failed"

# optional feature: autoheader
(autoheader --version) < /dev/null > /dev/null 2>&1
if [ $? -eq 0 ]; then
  echo "Running autoheader"
  autoheader --force \
      || error_and_exit "autoheader failed"
fi

echo "Running automake -Wall --copy --add-missing"
automake  -Wall --copy --add-missing \
    || error_and_exit "automake failed"


echo "Running autoconf"
    autoconf --force \
        || error_and_exit "autoconf failed"

echo "Success."
exit 0
