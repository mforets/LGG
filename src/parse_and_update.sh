#!/bin/sh

sage --preparse $1.sage
echo "Sage preparse ... OK"
mv $1.sage.py $1.py
