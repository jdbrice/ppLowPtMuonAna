## Running the Low Pt Dimuon Analysis

on macfrank the original code is here:
```bash
/home/jdb12/nwork/dimuonAna/ppLowPtMuonAna/
```

1. set up your environment
  - install [(mini)conda](https://conda.io/docs/user-guide/install/macos.html)
  - Add the following to your ~/.bashrc file
  ```bash
  export JDB_LIB="/macstar/star1/jdb12/RooBarb/"
  ```
  - At this point you should also see a line from the conda install like:
  ```bash
  export PATH="/macstar/star1/jdb12/vendor/conda2/bin:$PATH"
  ```
  but with a different path of course.
  - install scons with conda
  ```bash
  conda install scons
  ```
2. Copy and Build the code
  - copy the code to a new directory 
  - run the script (in bin/)
  ```bash
  ./fullbuild 
  ```
  to build the code, or just run 
  ```bash
  scons 
  ```
  from the root project directory
3. run analysis from inside bin:
```
./pairAna.app config/<config-file.xml>
```
skim.xml - reads the PicoDsts and produce histograms etc.
fit.xml - produce fits of the muon and pion peaks
pairs.xml - Do PID and make the invariant mass histograms based on the fits 


