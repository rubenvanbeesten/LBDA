# How to run the code

1. ...
2. ...

3. Run the sslp file data/[filename] by running the line
./lbdaexec data/sslp/[filename] -l -1000 -L 1.0 -B 0.0 -t 1 -T 3600 -m decom
Here,
* -L [l1 l2 ...]:   sequence of lambdas (indicating percentage share of each CVaR)
* -B [b1 b2 ...]:   ...
* -m [method]:      solution method
* -c [cut]:         type of cuts used (default: loose Benders cuts)
* -l [lb]:          lower bound of value [lb]
* -t [time]:           time limit of value t for each total run