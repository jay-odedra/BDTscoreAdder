# BDTScore Adder

install fastforest xgboost evaluator
```
git clone git@github.com:guitargeek/FastForest.git
cd FastForest
mkdir build
cd build
cmake ..
make
cp -P src/libfastforest.so* .
```
add bdtscoreadder repo

```
git clone git@github.com:jay-odedra/BDTscoreAdder.git
cd BDTscoreAdder
g++ -fPIC -std=c++11 Events.cc -o addbdtscore.exe -lfastforest -I ../FastForest/include/ -L ../FastForest/build/ `root-config --glibs --cflags`
./addbdtscore.exe inputfiles.txt outputfilename outputfiledir
```