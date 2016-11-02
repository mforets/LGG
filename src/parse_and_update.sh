#!/bin/sh

printf "sage: preparsing $1.sage ... "
sage --preparse $1.sage
printf "done\n"
mv $1.sage.py $1.py
