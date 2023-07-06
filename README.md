# BDTScore Adder
Install CMSSW
```
cmsrel CMSSW_13_1_0
cd CMSSW_13_1_0/src
cmsenv
```
install fastforest xgboost evaluator
```
git clone git@github.com:guitargeek/FastForest.git
```
build it
```
cd FastForest
mkdir build
cd build
cmake ..
make
cp -P src/libfastforest.so* .
```
add bdtscoreadder repo

```
cd ${CMSSW_BASE}/src
git clone git@github.com:jay-odedra/BDTscoreAdder.git
cd BDTscoreAdder
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CMSSW_BASE}/src/FastForest/build/
g++ -fPIC -std=c++11 Events.cc -o addbdtscore.exe -lfastforest -I ../FastForest/include/ -L ../FastForest/build/ `root-config --glibs --cflags`
./addbdtscore.exe inputfiles.txt outputfilename outputfiledir
```
XGBOOST MODELS SHOULD BE TRAINED USING OBJECTIVE = objective='binary:logitraw' AND SAVED USING booster.dump_model("model.txt")