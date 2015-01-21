#!/bin/sh

rm -rf autom4te.cache
rm -f aclocal.m4 adm_local/ltmain.sh

echo "Running aclocal..."    ;
aclocal --force -I adm_local || exit 1
echo "Running autoheader..." ; autoheader --force -I adm_local            || exit 1
echo "Running autoconf..."   ; autoconf --force                    || exit 1
echo "Running libtoolize..." ; libtoolize --copy --force           || exit 1
echo "Running automake..."   ; automake --add-missing --copy       || exit 1
