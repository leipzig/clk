#!/bin/bash
export HTSLIB_VERSION='1.3.2'
#for zlib
export CFLAGS="-I$PREFIX/include"
export LDFLAGS="-L$PREFIX/lib"

mkdir -p $PREFIX/bin
echo "#!/usr/bin/env python3" > $PREFIX/bin/rMATS-ISO.py
cat $SRC_DIR/rMATS-ISO.py >> $PREFIX/bin/rMATS-ISO.py
cp utils.py $PREFIX/bin
chmod +x $PREFIX/bin/rMATS-ISO.py

mkdir -p $PREFIX/config
cp lr2rmats/config.template.yaml lr2rmats/Snakefile $PREFIX/config/

cd IsoModule
# Remove rdynamic which can cause build issues on OSX
# https://sourceforge.net/p/samtools/mailman/message/34699333/
sed -i.bak 's/ -rdynamic//g' Makefile
sed -i.bak 's/ -rdynamic//g' htslib-${HTSLIB_VERSION}/configure

export CPPFLAGS="-I$PREFIX/include"
export LDFLAGS="-L$PREFIX/lib"

# Problem with ncurses from default channel we now get in bioconda so skip tview
# https://github.com/samtools/samtools/issues/577
cd htslib-${HTSLIB_VERSION}
#if [ "$(uname)" == "Darwin" ]; then
  # there is some problem on High Sierra where inttypes.h is not found, leading to problems compiling lzma.h
#  ./configure --prefix=$PREFIX --enable-libcurl --disable-lzma CFLAGS="-I$PREFIX/include" LDFLAGS="-L$PREFIX/lib" 
#else
# --enable-libcurl
./configure --prefix=$PREFIX CFLAGS="-I$PREFIX/include" LDFLAGS="-L$PREFIX/lib"
#fi
make
make install

#to IsoModule
cd ..
ln -s htslib-${HTSLIB_VERSION} htslib
make CFLAGS='-Wall -O2 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -I${PREFIX}/include/' LDFLAGS="-L$PREFIX/lib"
cp IsoModule $PREFIX/bin
#to top level


cd $SRC_DIR

cd lr2rmats
ln -s ../IsoModule/htslib .
make CFLAGS='-Wall -O2 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -I${PREFIX}/include/' LDFLAGS="-L$PREFIX/lib"
cp bin/lr2rmats $PREFIX/bin

cd $SRC_DIR
make all

echo "done installing"


# WARNING (rmats-iso,bin/IsoModule): did not find - or even know where to look for: /lib64/libc.so.6
#   INFO (rmats-iso,bin/IsoModule): Needed DSO lib/libz.so.1 found in conda-forge::zlib-1.2.11-h470a237_3
# WARNING (rmats-iso,bin/IsoModule): did not find - or even know where to look for: /lib64/libm.so.6
# WARNING (rmats-iso,bin/IsoModule): did not find - or even know where to look for: /lib64/libpthread.so.0
# WARNING (rmats-iso,bin/bgzip): did not find - or even know where to look for: /lib64/libc.so.6
#   INFO (rmats-iso,bin/bgzip): Needed DSO lib/libz.so.1 found in conda-forge::zlib-1.2.11-h470a237_3
# WARNING (rmats-iso,bin/bgzip): did not find - or even know where to look for: /lib64/libpthread.so.0
# WARNING (rmats-iso,bin/tabix): did not find - or even know where to look for: /lib64/libc.so.6
#   INFO (rmats-iso,bin/tabix): Needed DSO lib/libz.so.1 found in conda-forge::zlib-1.2.11-h470a237_3
# WARNING (rmats-iso,bin/tabix): did not find - or even know where to look for: /lib64/libpthread.so.0
#   INFO (rmats-iso,lib/libhts.so.1.3.2): Needed DSO lib/libz.so.1 found in conda-forge::zlib-1.2.11-h470a237_3
# WARNING (rmats-iso,lib/libhts.so.1.3.2): did not find - or even know where to look for: /lib64/libc.so.6
# WARNING (rmats-iso,lib/libhts.so.1.3.2): did not find - or even know where to look for: /lib64/libm.so.6
# WARNING (rmats-iso,lib/libhts.so.1.3.2): did not find - or even know where to look for: /lib64/libpthread.so.0
# WARNING (rmats-iso,bin/htsfile): did not find - or even know where to look for: /lib64/libc.so.6
#   INFO (rmats-iso,bin/htsfile): Needed DSO lib/libz.so.1 found in conda-forge::zlib-1.2.11-h470a237_3
# WARNING (rmats-iso,bin/htsfile): did not find - or even know where to look for: /lib64/libpthread.so.0