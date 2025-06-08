#!/bin/bash 
# used to update Makefile with the last git version available
MYVERSION="distribution-version"
if [ -x "$(command -v git)" ]; then
    MYVERSION=`git describe --abbrev=7 --dirty --always --tag --long 2> /dev/null`
    if [ -z "$MYVERSION" ] ; then
	MYVERSION="Distribution-version"
    fi
fi

echo $MYVERSION
