#!/bin/sh

printf "sage: preparsing lgg.sage ... "
sage --preparse lgg.sage
printf "done\n"
mv lgg.sage.py lgg.py
