### Minimum requirements to compile:
C++ compiler that supports C++11 standard


### Optional libraries for better performance:
BLAS, LAPACK, FFTW

### Instructions:

```bash
git clone https://github.com/dmalhotra/test-slenderbody.git
cd test-slenderbody
git submodule init
git submodule update
# update Makefile as necessary
make
./bin/slenderbody
```

For tangle test do:

```bash
make && ./bin/test-unitDLP-tangle && paraview tangleUerr.pvtu
```
