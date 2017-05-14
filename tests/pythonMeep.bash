#/bin/bash

#
# test for successful loading of python meep module
#
cd ../python
ln -f -s .libs/_meep.so
python -c 'import meep'
