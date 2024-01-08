# PASS, Presorted Alignment for Short Sequences

The tool we used to build Affymetrix custom CDF for two decades.

### Install
```
git clone https://github.com/daimh/pass.git
cd pass
autoreconf -i
./configure
make
make install # or copy the binary to one of your PATHS, such as 'cp src/pass ~/bin'
pass --help
```

## Acknowledgments

Fan Meng, Ph.D., Research Associate Prof, Psychiatry, UMICH

David Philip Kreil, PH.D., Assoc. Prof, Universität für Bodenkultur Wien

Huda Akil, Ph.D., Director of MNI, UMICH

Stanley J. Watson, M.D., Ph.D., Director of MNI, UMICH


### License

MIT License

Copyright © 2002-2024 University of Michigan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
